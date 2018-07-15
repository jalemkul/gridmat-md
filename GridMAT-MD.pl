#!/usr/bin/perl

#########################################################################
#       	This is GridMAT-MD, version 2.0		       		            #
#		Release date: 04/18/2013				                        #
#									                                    #
#		For bug fixes and release notes				                    #
#		please see the following site:				                    #
#	http://www.bevanlab.biochem.vt.edu/GridMAT-MD/bugs.html		        #
#									                                    #
#########################################################################



##### Program Begin #####
use strict;
my $start_time = time();



##### Read parameter file #####
unless (@ARGV){ die "Usage: ./$0 <param file>\n"; }

my $param_file = @ARGV[0];
open (INPUT, $param_file) or die "It seems \"$param_file\" does not exist!\n";
my @param_file = (<INPUT>);
close (INPUT);

my @param_file_short = grep !/^#/, @param_file;
my %param;
foreach (@param_file_short){
	#my ($key, $value) = split(/\s+/, $_);
	my ($key, $value) = split(' ', $_);
	$param{$key} = $value;
}

print "\nReading from \"$param_file\"...\n";
print "You defined the coordinate file as $param{coord_file} \n";
print "You specified that this file contains $param{num_frames} frame(s) \n";

for (my $i=1; $i<$param{num_lipid_types}+1; $i++){
	my $res = "resname".$i;
	my $atom = "atomname".$i;
	print "You defined a lipid residue name as $param{$res} (atom(s): $param{$atom}) \n";
}

print "You defined the solvent as $param{solvent} \n";
print "You defined the ions as $param{ions} \n";



##### Calculate bilayer thickness and areas; write output files #####
if ($param{num_frames} > 1){

	for (my $i=0; $i<$param{num_frames}; $i++){
		print "\nRunning GridMAT on frame ".($i+1).":\n";
		&gridmat($i);
	}

} else {
	&gridmat;
}



##### Program End #####
print "\nThank you for using this program.  Please read and cite the following reference:\n";
print "\nW. J. Allen, J. A. Lemkul, and D. R. Bevan. (2009) \"GridMAT-MD: A Grid-based Membrane\n";
print "Analysis Tool for Use With Molecular Dynamics.\" J. Comput. Chem., 30 (12), 1952-1958\n";

my $end_time = time();
print "\nTime taken was: ".($end_time-$start_time)." seconds.\n";
exit();



#########################################################################
#########################################################################
## Subroutines


sub gridmat{

	# Read coordinate file
	my @bilayer;
	my $num_lines = `wc -l $param{coord_file}`;
	my $num_lines_each = $num_lines / $param{num_frames};
	my $start = 1 + $_[0]*$num_lines_each;
	my $end = $start + $num_lines_each - 1;

	open (INPUT2, $param{coord_file}) or die "It seems \"$param{coord_file}\" does not exist!\n";
	while (<INPUT2>){
		if ($. > $end){ last; }
		if ($. >= $start){ push (@bilayer, $_); }
	}
	close (INPUT2);



	# Store atom information in appropriate format (@bilayer_fixed)
	my (@bilayer_fixed, $box_vectors);

	if ($param{file_type} =~ /pdb/){

		my @atoms = grep /^ATOM|^HETATM/, @bilayer;
		my @cryst = grep /^CRYST1/, @bilayer;
		my ($temp1, $temp2, $temp3, $temp4, $temp5) = split(' ', $cryst[0]);
		$box_vectors = ($temp2/10)."\t".($temp3/10)."\t".($temp4/10)."\n";

		foreach (@atoms){

			my @line = split(//, $_);

			my $z_coor = sprintf("%.3f", (join("", splice(@line, 46, 8)))/10);
			my $y_coor = sprintf("%.3f", (join("", splice(@line, 38, 8)))/10);
			my $x_coor = sprintf("%.3f", (join("", splice(@line, 30, 8)))/10);
			my $r_num  = join("", splice(@line, 22, 4));
			my $r_name = join("", splice(@line, 17, 4));
			my $a_name = join("", splice(@line, 12, 4));
			my $a_num  = join("", splice(@line, 6, 5));

			my $new_line = ("$r_num $r_name  $a_name  $a_num  $x_coor  $y_coor  $z_coor");
			push(@bilayer_fixed, $new_line."\n");
		}

	} elsif ($param{file_type} =~ /gro/){

		shift @bilayer;
		shift @bilayer;
		$box_vectors = pop @bilayer;

		foreach (@bilayer){
			my $newline = pack("A6 A5 A7 A29", unpack("A5 A4 A6 A29", $_));
			push (@bilayer_fixed, $newline."\n");
		}

	} else {

		print "\nParameter file_type=$param{file_type} not recognized! Check your param file and coordinate file!\n";
		exit();
	}



	# Iterate through the coordinate file and pull out the atoms that match any input resnameX / atomnameX
	my @pre_reference;
	my $greatest_reference = 0;


	for (my $i=1; $i<$param{num_lipid_types}+1; $i++){
		my $res = "resname".$i;
		my $pre_atom = "atomname".$i;
		my @atomnames = split(",", $param{$pre_atom});


		foreach (@bilayer_fixed){
			my ($b_resid, $b_resname, $b_atomname, $b_atomid, $bx, $by, $bz) = split(' ', $_);

			if ($b_resname eq($param{$res})){
				foreach my $atom (@atomnames){
					if ($atom eq($b_atomname)){
						push (@pre_reference, $b_resid." ".$b_atomname." ".$bx." ".$by." ".$bz."\n");
						if ($b_resid > $greatest_reference){
							$greatest_reference = $b_resid;
						}
					}
				}
			}
		}
	}



	# Find the center-of-mass of reference atoms, if multiple are declared
	my @reference_atoms;
	for (my $i=0; $i<($greatest_reference+1); $i++){
		my ($x_tot, $y_tot, $z_tot, $x_ave, $y_ave, $z_ave, $flag) = 0;
	
		foreach (@pre_reference){
			my ($b_resid, $b_atomname, $bx, $by, $bz) = split(' ', $_);
	
			if ($b_resid == $i){
				$x_tot += $bx;
				$y_tot += $by;
				$z_tot += $bz;
				$flag += 1;
			}
		}
	
		if ($flag > 0){
			$x_ave = sprintf("%.3f", ($x_tot / $flag));
			$y_ave = sprintf("%.3f", ($y_tot / $flag));
			$z_ave = sprintf("%.3f", ($z_tot / $flag));
	
			push(@reference_atoms, $i." ".$x_ave." ".$y_ave." ".$z_ave."\n");
		}
	}



	# Figure out box size
	my ($ref_z_max, $ref_z_min) = &limits(\@reference_atoms, 4);
	my $z_middle = ($ref_z_max + $ref_z_min)/2;
	my ($ref_x_max, $ref_x_min, $ref_y_max, $ref_y_min);
	my ($transfer_x, $transfer_y);
	my $box_size_source = "a line in the coord file";

	if ($param{box_size} eq("vectors")){

		if (exists($param{override_vectors})){
			my ($vec1, $vec2, $vec3) = split(",", $param{override_vectors});
			$box_vectors = ($vec1." ".$vec2." ".$vec3);
			$box_size_source = "the parameter file";
		}


		($transfer_x, $transfer_y) = split(' ', $box_vectors);
	
		($ref_x_max) = &limits(\@reference_atoms, 2);
		($ref_y_max) = &limits(\@reference_atoms, 3);
		$ref_x_min = $ref_x_max - $transfer_x;
		$ref_y_min = $ref_y_max - $transfer_y;

	} elsif ($param{box_size} eq("solvent")){
		my @solvent_atoms;	
	
		foreach (@bilayer_fixed){
			my ($resid, $resname, $atomname, $atomid, $x, $y, $z) = split(' ', $_);
			if ($resname eq($param{solvent})){
				push(@solvent_atoms, $_);
			}
		}
	
		($ref_x_max, $ref_x_min) = &limits(\@solvent_atoms, 5);
		($ref_y_max, $ref_y_min) = &limits(\@solvent_atoms, 6);
	
		$transfer_x = $ref_x_max - $ref_x_min;
		$transfer_y = $ref_y_max - $ref_y_min;	
		$box_size_source = "the solvent coordinates";

	} else {

		print "\nParameter box_size = $param{box_size} not recognized! Check your param file.\n";
		exit();
	}	

	print "\nDeconstructing lipid bilayer...\n";
	print "Lower X limit: $ref_x_min   Upper X limit: $ref_x_max\n";
	print "Lower Y limit: $ref_y_min   Upper Y limit: $ref_y_max\n";
	print "Cross sectional area (box size) was determined from: $box_size_source\n";
	print "Cross sectional area of the system: $transfer_x x $transfer_y nanometers\n";
	print "Lower Z limit: $ref_z_min   Upper Z limit: $ref_z_max\n";
	print "The middle (in the Z-direction) is $z_middle\n";



	# Divide reference points into a top and bottom leaflet
	my (@top_leaflet_ref, @bottom_leaflet_ref);
	foreach (@reference_atoms){
		my ($resid, $refx, $refy, $refz) = split(' ', $_);
		if ($refz > $z_middle){
			push(@top_leaflet_ref, $_);
		} else {
			push(@bottom_leaflet_ref, $_);
		}
	}
	
	my ($top_z_max, $top_z_min) = &limits(\@top_leaflet_ref, 4);
	my ($bottom_z_max, $bottom_z_min) = &limits(\@bottom_leaflet_ref, 4);
	
	print "In the top leaflet, the Z values range from $top_z_min to $top_z_max\n";
	print "In the bottom leaflet, the Z values range from $bottom_z_min to $bottom_z_max\n";



	# Set up periodic images of the bilayer
	print "\nSimulating periodic boundary conditions...";
	my @pbc_atoms;
	my (@a_1A, @a_2A, @a_3A, @a_1B, @a_3B, @a_1C, @a_2C, @a_3C);
	
	foreach (@reference_atoms){
		my ($resid, $oldx, $oldy, $oldz) = split(' ', $_);
	
		my $x_A = sprintf("%.3f", ($oldx-$transfer_x));
		my $x_B = $oldx;
		my $x_C = sprintf("%.3f", ($oldx+$transfer_x));
	
		my $y_1 = sprintf("%.3f", ($oldy+$transfer_y));
		my $y_2 = $oldy;
		my $y_3 = sprintf("%.3f", ($oldy-$transfer_y));
	
		push(@a_1A, $resid." ".$x_A." ".$y_1." ".$oldz."\n");
		push(@a_2A, $resid." ".$x_A." ".$y_2." ".$oldz."\n");
		push(@a_3A, $resid." ".$x_A." ".$y_3." ".$oldz."\n");
	
		push(@a_1B, $resid." ".$x_B." ".$y_1." ".$oldz."\n");
		push(@a_3B, $resid." ".$x_B." ".$y_3." ".$oldz."\n");
	
		push(@a_1C, $resid." ".$x_C." ".$y_1." ".$oldz."\n");
		push(@a_2C, $resid." ".$x_C." ".$y_2." ".$oldz."\n");
		push(@a_3C, $resid." ".$x_C." ".$y_3." ".$oldz."\n");
	}
	
	@pbc_atoms = (@a_1A, @a_1B, @a_1C, @a_2A, @reference_atoms, @a_2C, @a_3A, @a_3B, @a_3C);
	print "Done\nDividing the periodic array into a top and bottom leaflet...";
	
	my (@top_leaflet_ref_pbc, @bottom_leaflet_ref_pbc);
	foreach (@pbc_atoms){
		my ($resid, $bx, $by, $bz) = split(' ', $_);
		if ($bz > $z_middle){
			push(@top_leaflet_ref_pbc, $_);
		} else {
			push(@bottom_leaflet_ref_pbc, $_);
		}
	}
	
	print "Done\n";



	# Look for protein atoms that fall within the head groups of the lipids
	my (@bilayer_protein, @top_offenders, @top_offenders_final, @bottom_offenders, @bottom_offenders_final);
	
	if ($param{protein} eq("yes")){
		print "\nLooking for offending protein atoms...\n";

		for (my $i=1; $i<$param{num_lipid_types}+1; $i++){
			my $res = "resname".$i;

			if ($i == 1){
				@bilayer_protein = grep(!/$param{$res}/, @bilayer_fixed);
			} else {
				@bilayer_protein = grep(!/$param{$res}/, @bilayer_protein);
			}
		}
	
		if (exists $param{solvent}){
			@bilayer_protein = grep(!/$param{solvent}/, @bilayer_protein);
		}
	
		if (exists $param{ions}){
			my @ions = split(",", $param{ions});
			foreach my $ion (@ions){
				@bilayer_protein = grep(!/$ion/, @bilayer_protein);
			}
		}

		foreach (@bilayer_protein){
			my ($resid, $resname, $atomname, $atomid, $px, $py, $pz) = split(' ', $_);
			if ($pz <= $top_z_max && $pz >= $top_z_min){
				push (@top_offenders, "protein ".$px." ".$py." ".$pz."\n");
			}
			if ($pz <= $bottom_z_max && $pz >= $bottom_z_min){
				push (@bottom_offenders, "protein ".$px." ".$py." ".$pz."\n");
			}
		}
	
		foreach (@top_offenders){
			my ($resid, $px, $py, $pz) = split(' ', $_);
			my @top_neighbors = "";
	
			foreach (@top_leaflet_ref_pbc){
				my ($ref_resid, $ref_x, $ref_y, $ref_z) = split(' ', $_);
	
				if (sqrt(($px - $ref_x)**2 + ($py - $ref_y)**2) <= $param{precision}){
					push(@top_neighbors, $ref_z);
				}
			}
			my ($z_greatest, $z_lowest) = &limits(\@top_neighbors, 1);
	
			if ($pz >= $z_lowest && $pz <= $z_greatest){
				push (@top_offenders_final, "protein ".$px." ".$py." ".$pz."\n");
			}
		}
	
		foreach (@bottom_offenders){
			my ($resid, $px, $py, $pz) = split(' ', $_);
			my @bottom_neighbors = "";
	
			foreach (@bottom_leaflet_ref_pbc){
				my ($ref_resid, $ref_x, $ref_y, $ref_z) = split(' ', $_);
	
				if (sqrt(($px - $ref_x)**2 + ($py - $ref_y)**2) <= $param{precision}){
					push(@bottom_neighbors, $ref_z);
				}
			}
			my ($z_greatest, $z_lowest) = &limits(\@bottom_neighbors, 1);
	
			if ($pz >= $z_lowest && $pz <= $z_greatest){
				push (@bottom_offenders_final, "protein ".$px." ".$py." ".$pz."\n");
			}
		}
	
		print "There are ".scalar(@top_offenders_final)." protein atoms within the headgroups of the top leaflet\n";
		print "There are ".scalar(@bottom_offenders_final)." protein atoms within the headgroups of the bottom leaflet\n";
	}



	# Set up periodic images of the bilayer + protein
	my @reference_atoms_protein = (@reference_atoms, @top_offenders_final, @bottom_offenders_final);
	my (@pbc_atoms_protein, @top_leaflet_protein_pbc, @bottom_leaflet_protein_pbc);
	my (@p_1A, @p_2A, @p_3A, @p_1B, @p_3B, @p_1C, @p_2C, @p_3C);
	
	if ($param{protein} eq("yes")){
	
		print "\nSimulating periodic boundary conditions for the protein atoms...";
	
		foreach (@reference_atoms_protein){
			my ($resid, $oldx, $oldy, $oldz) = split(' ', $_);
	
			my $x_A = sprintf("%.3f", ($oldx-$transfer_x));
			my $x_B = $oldx;
			my $x_C = sprintf("%.3f", ($oldx+$transfer_x));
	
			my $y_1 = sprintf("%.3f", ($oldy+$transfer_y));
			my $y_2 = $oldy;
			my $y_3 = sprintf("%.3f", ($oldy-$transfer_y));
	
			push(@p_1A, $resid." ".$x_A." ".$y_1." ".$oldz."\n");
			push(@p_2A, $resid." ".$x_A." ".$y_2." ".$oldz."\n");
			push(@p_3A, $resid." ".$x_A." ".$y_3." ".$oldz."\n");
	
			push(@p_1B, $resid." ".$x_B." ".$y_1." ".$oldz."\n");
			push(@p_3B, $resid." ".$x_B." ".$y_3." ".$oldz."\n");
	
			push(@p_1C, $resid." ".$x_C." ".$y_1." ".$oldz."\n");
			push(@p_2C, $resid." ".$x_C." ".$y_2." ".$oldz."\n");
			push(@p_3C, $resid." ".$x_C." ".$y_3." ".$oldz."\n");
		}
	
		@pbc_atoms_protein = (@p_1A, @p_1B, @p_1C, @p_2A, @reference_atoms_protein, @p_2C, @p_3A, @p_3B, @p_3C);
		print "Done\nDividing the periodic array into a top and bottom leaflet...";
	
		foreach (@pbc_atoms_protein){
			my ($resid, $bx, $by, $bz) = split(' ', $_);
			if ($bz > $z_middle){
				push(@top_leaflet_protein_pbc, $_);
			} else {
				push(@bottom_leaflet_protein_pbc, $_);
			}
		}
	
		print "Done\n";
	}



	# Generate the grid
	print "\nGenerating the grid...\n";
	my (@grid_final, @grid_x, @grid_y, $grid_x, $grid_y);
	
	if ($param{conserve_ratio} eq("yes")){
	
		if ($transfer_x >= $transfer_y){
			$grid_x = $param{grid};
			my $interval_x = ($transfer_x)/($grid_x-1);
			$grid_y = sprintf("%.0f", ($transfer_y/$interval_x) + 1);
			my $interval_y = ($transfer_y)/($grid_y-1);
	
			print "Your system is bigger in the X-direction\n";
			print "There are ".$grid_x." grid points in the X direction, spaced every ";
			printf ("%.5f", $interval_x);
			print " nanometers\n";
			print "There are ".$grid_y." grid points in the Y direction, spaced every ";
			printf ("%.5f", $interval_y);
			print " nanometers\n";
			print "Note: the intervals may not be exactly the same in order to have a whole number of grid points\n";
	
			for (my $n=0; $n<=($grid_x-1); $n++){
				my $grid_point = sprintf("%.3f", $ref_x_min + ($interval_x * $n));
				push(@grid_x, $grid_point."\n");
			}
	
			for (my $n=0; $n<=($grid_y-1); $n++){
				my $grid_point = sprintf("%.3f", $ref_y_min + ($interval_y * $n));
				push(@grid_y, $grid_point."\n");
			}
	
			@grid_final = &makegrid(\@grid_x, \@grid_y);
		}
	
		if ($transfer_y > $transfer_x){
			$grid_y = $param{grid};
			my $interval_y = ($transfer_y)/($grid_y-1);
			$grid_x = sprintf("%.0f", ($transfer_x/$interval_y) + 1);
			my $interval_x = ($transfer_x)/($grid_x-1);
	
			print "Your system is bigger in the Y-direction\n";
			print "There are ".$grid_x." grid points in the X direction, spaced every ";
			printf ("%.5f", $interval_x);
			print " nanometers\n";
			print "There are ".$grid_y." grid points in the Y direction, spaced every ";
			printf ("%.5f", $interval_y);
			print " nanometers\n";
			print "Note: the intervals may not be exactly the same in order to have a whole number of grid points\n";
	
			for (my $n=0; $n<=($grid_x-1); $n++){
				my $grid_point = sprintf("%.3f", $ref_x_min + ($interval_x * $n));
				push(@grid_x, $grid_point."\n");
			}
	
			for (my $n=0; $n<=($grid_y-1); $n++){
				my $grid_point = sprintf("%.3f", $ref_y_min + ($interval_y * $n));
				push(@grid_y, $grid_point."\n");
			}
	
			@grid_final = &makegrid(\@grid_x, \@grid_y);
		}
	
	} else {
	
		$grid_x = $param{grid};
		$grid_y = $param{grid};
		my $interval_x = ($transfer_x)/($grid_x-1);
		my $interval_y = ($transfer_y)/($grid_y-1);
	
		print "You defined the grid as ".$grid_x."x".$grid_y." points\n";
		printf ("X-grid points are spaced every %.5f nanometers\n", $interval_x);
		printf ("Y-grid points are spaced every %.5f nanometers\n", $interval_y);
	
		for (my $n=0; $n<=($grid_x-1); $n++){
			my $grid_point = sprintf("%.3f", $ref_x_min + ($interval_x * $n));
			push(@grid_x, $grid_point."\n");
		}
	
		for (my $n=0; $n<=($grid_y-1); $n++){
			my $grid_point = sprintf("%.3f", $ref_y_min + ($interval_y * $n));
			push(@grid_y, $grid_point."\n");
		}
	
		@grid_final = &makegrid(\@grid_x, \@grid_y);
	}



	# Copy the correct coordinates to the top_leaflet / bottom_leaflet arrays
	print "\nAnalyzing the bilayer...\n";
	my (@top_leaflet, @bottom_leaflet);
	
	if ($param{protein} eq("yes")){
		@top_leaflet = @top_leaflet_protein_pbc;
		@bottom_leaflet = @bottom_leaflet_protein_pbc;
	} else {
		@top_leaflet = @top_leaflet_ref_pbc;
		@bottom_leaflet = @bottom_leaflet_ref_pbc;
	}



	# Prepare the output filename prefix
	my $output_filename;
	$output_filename .= $param{output_prefix};
	$output_filename .= ".frame";
	$output_filename .= ($_[0]+1);
	$output_filename .= ".";
	$output_filename .= $grid_x."x".$grid_y;
	


	# Calculate the thicknesses from the perspective of the top leaflet
	my @top_leaflet_z_values;
	
	{
		my $closest = 1000000;
		my $closest_z;
	
		foreach my $m (@top_leaflet){
			my ($resid_top, $topx, $topy, $topz) = split(' ', $m);
	
			foreach my $n (@bottom_leaflet){
				my ($resid_bottom, $bottomx, $bottomy, $bottomz) = split(' ', $n);
				my $dist = (sqrt(($topx-$bottomx)**2 + ($topy-$bottomy)**2));
	
				if ($dist < $closest){
					$closest = $dist;
					$closest_z = ($topz-$bottomz);
				}
			}
			push (@top_leaflet_z_values, $resid_top." ".$topx." ".$topy." ".$closest_z."\n");
			$closest = 1000000;
		}
	}
	
	my (@top_final_a, @top_final_z);
	
	{
		my $closest = 1000000;
		my ($closest_z, $lowest_resid);
	
		foreach my $m (@grid_final){
			my ($gridx, $gridy) = split(' ', $m);
	
			foreach my $n (@top_leaflet_z_values){
				my ($resid, $topx, $topy, $topz) = split(' ', $n);
				my $dist = (sqrt(($gridx-$topx)**2 + ($gridy-$topy)**2));
	
				if ($dist < $closest){
					$closest = $dist;
					$closest_z = $topz;
					$lowest_resid = $resid;
				}
			}
			
			if ($lowest_resid eq("protein")){
				$closest_z = $param{P_value};
			}
	
			push(@top_final_a, $lowest_resid." ".$closest_z."\n");
			push(@top_final_z, $closest_z."\n");
			$closest = 1000000;
		}
	}
	
	# Write the top leaflet thicknesses to file
	if ($param{thickness} eq("yes")){

		if ($param{output_format} eq("column")){

			print "The top leaflet thickness will be printed to ".$output_filename.".top_thickness.dat\n";
			open(OUTFILE1, ">", $output_filename.".top_thickness.dat");
			print OUTFILE1 @top_final_z;
			close(OUTFILE1);

		} elsif ($param{output_format} eq("matrix")){

			&convert_to_matrix(\@top_final_z, $grid_x, $grid_y, $output_filename, ".top_thickness");

		} elsif ($param{output_format} eq("vector")){

			&convert_to_vector(\@top_final_z, $grid_x, $grid_y, $output_filename, ".top_thickness");

		} else {

			print "\nParameter output_format = $param{output_format} not recognized! Check your param file.\n";
			exit();
		}
	}
	


	# Calculate the thicknesses from the perspective of the bottom leaflet
	my @bottom_leaflet_z_values;
	
	{
		my $closest = 1000000;
		my $closest_z;
	
		foreach my $m (@bottom_leaflet){
			my ($resid_bottom, $bottomx, $bottomy, $bottomz) = split(' ', $m);
	
			foreach my $n (@top_leaflet){
				my ($resid_top, $topx, $topy, $topz) = split(' ', $n);
				my $dist = (sqrt(($bottomx-$topx)**2 + ($bottomy-$topy)**2));
	
				if ($dist < $closest){
					$closest = $dist;
					$closest_z = ($topz-$bottomz);
				}
			}
			push (@bottom_leaflet_z_values, $resid_bottom." ".$bottomx." ".$bottomy." ".$closest_z."\n");
			$closest = 1000000;
		}
	}
	
	my (@bottom_final_a, @bottom_final_z);
	
	{
		my $closest = 1000000;
		my ($closest_z, $lowest_resid);
	
		foreach my $m (@grid_final){
			my ($gridx, $gridy) = split(' ', $m);
	
			foreach my $n (@bottom_leaflet_z_values){
				my ($resid, $bottomx, $bottomy, $bottomz) = split(' ', $n);
				my $dist = (sqrt(($gridx-$bottomx)**2 + ($gridy-$bottomy)**2));
	
				if ($dist < $closest){
					$closest = $dist;
					$closest_z = $bottomz;
					$lowest_resid = $resid;
				}
			}
			
			if ($lowest_resid eq("protein")){
				$closest_z = $param{P_value};
			}
	
			push(@bottom_final_a, $lowest_resid." ".$closest_z."\n");
			push(@bottom_final_z, $closest_z."\n");
			$closest = 1000000;
		}
	}
	
	# Write the bottom leaflet thicknesses to file
	if ($param{thickness} eq("yes")){

		if ($param{output_format} eq("column")){

			print "The bottom leaflet thickness will be printed to ".$output_filename.".bottom_thickness.dat\n";
			open(OUTFILE2, ">", $output_filename.".bottom_thickness.dat");
			print OUTFILE2 @bottom_final_z;
			close(OUTFILE2);

		} elsif ($param{output_format} eq("matrix")){

			&convert_to_matrix(\@top_final_z, $grid_x, $grid_y, $output_filename, ".bottom_thickness");

		} elsif ($param{output_format} eq("vector")){

			&convert_to_vector(\@top_final_z, $grid_x, $grid_y, $output_filename, ".bottom_thickness");

		}
	}



	# Calculate average thicknesses
	my @average_final;
	
	{
		my $average_value;
		for (my $n=0; $n<(scalar(@grid_final)); $n++){
			$average_value = ($top_final_z["$n"] + $bottom_final_z["$n"])/2;
			push (@average_final, $average_value."\n");
		}
	}
	
	# Write the average thicknesses to file
	if ($param{thickness} eq("yes")){

		if ($param{output_format} eq("column")){

			print "The average bilayer thickness will be printed to ".$output_filename.".average_thickness.dat\n";
			open(OUTFILE2, ">", $output_filename.".average_thickness.dat");
			print OUTFILE2 @average_final;
			close(OUTFILE2);

		} elsif ($param{output_format} eq("matrix")){

			&convert_to_matrix(\@top_final_z, $grid_x, $grid_y, $output_filename, ".average_thickness");

		} elsif ($param{output_format} eq("vector")){

			&convert_to_vector(\@top_final_z, $grid_x, $grid_y, $output_filename, ".average_thickness");

		}

		if ($param{protein} eq("yes")){
			print "Note: the average thickness file does not make much sense if protein = yes. See the manual.\n";
		}
	}



	# Print APLs just based on number of lipids and cross-sectional surface area
	print "\nCalculating area per lipid head group...\n";
	
	my $area = $transfer_x*$transfer_y*100;
	my $area_per_lipid_top = $area/(scalar(@top_leaflet_ref));
	my $area_per_lipid_bottom = $area/(scalar(@bottom_leaflet_ref));
	
	printf ("The lateral area of the system is %.5f sq. Angstroms (per side)\n", $area);
	print "When you don't account for any protein atoms:\n";
	printf ("    The average area per lipid in the top leaflet is %.5f sq. Angstroms\n", $area_per_lipid_top);
	printf ("    The average area per lipid in the bottom leaflet is %.5f sq. Angstroms\n", $area_per_lipid_bottom);
	
	if ($param{protein} eq("yes")){
		print "When you do take the protein atoms into account:\n";
	}
	


	# Calculate specific top leaflet areas
	my (@top_areas_sorted, $new_area_per_lipid_top);
	
	{
		my (%count_top, @count_top, @top_areas, $top_area_protein);
		push(@top_areas, "Resid\tZ-value\tArea (sq. Ang.)\n");
		chomp(@top_final_a);
	
		foreach (@top_final_a){
			if (exists $count_top{"$_"}){
				$count_top{"$_"}++;
			} else {
				$count_top{"$_"} = 1;
			}
		}
		while (my ($key, $value) = each %count_top){
			push(@count_top, $key." ".$value."\n");
		}
	
		foreach (@count_top){
			my ($resid, $z_value, $count) = split(' ', $_);
			my $new_area = (($count/($grid_x*$grid_y))*$area);
			push(@top_areas, $resid."\t".$z_value."\t".$new_area."\n");
			if ($resid eq("protein")){
				$top_area_protein = $new_area;
			}
		}
	
		my %sort_top = map {$_, 1} @top_areas;
		@top_areas_sorted = sort {$a <=> $b} (keys %sort_top);
	
		if ($param{protein} eq("yes")){
			my $new_area_top = $area - $top_area_protein;
			$new_area_per_lipid_top = $new_area_top/(scalar(@top_leaflet_ref));
			printf ("    The new area per lipid in the top leaflet is %.5f sq. Angstroms\n", $new_area_per_lipid_top);
			unshift(@top_areas_sorted, "Ave APL = $new_area_per_lipid_top sq. Ang\n\n");
		} else {
			unshift(@top_areas_sorted, "Ave APL = $area_per_lipid_top sq. Ang\n\n");
		}
	}



	# Calculate specific bottom leaflet areas
	my (@bottom_areas_sorted, $new_area_per_lipid_bottom);
	
	{
		my (%count_bottom, @count_bottom, @bottom_areas, $bottom_area_protein);
		push(@bottom_areas, "Resid\tZ-value\tArea (sq. Ang.)\n");
		chomp(@bottom_final_a);
	
		foreach (@bottom_final_a){
			if (exists $count_bottom{"$_"}){
				$count_bottom{"$_"}++;
			} else {
				$count_bottom{"$_"} = 1;
			}
		}
		while (my ($key, $value) = each %count_bottom){
			push(@count_bottom, $key." ".$value."\n");
		}
	
		foreach (@count_bottom){
			my ($resid, $z_value, $count) = split(' ', $_);
			my $new_area = (($count/($grid_x*$grid_y))*$area);
			push(@bottom_areas, $resid."\t".$z_value."\t".$new_area."\n");
			if ($resid eq("protein")){
				$bottom_area_protein = $new_area;
			}
		}
	
		my %sort_bottom = map {$_, 1} @bottom_areas;
		@bottom_areas_sorted = sort {$a <=> $b} (keys %sort_bottom);
	
		if ($param{protein} eq("yes")){
			my $new_area_bottom = $area - $bottom_area_protein;
			$new_area_per_lipid_bottom = $new_area_bottom/(scalar(@bottom_leaflet_ref));
			printf ("    The new area per lipid in the bottom leaflet is %.5f sq. Angstroms\n", $new_area_per_lipid_bottom);
			unshift(@bottom_areas_sorted, "Ave APL = $new_area_per_lipid_bottom sq. Ang\n\n");
		} else {
			unshift(@bottom_areas_sorted, "Ave APL = $area_per_lipid_bottom sq. Ang\n\n");
		}
	}



	# Print top leaflet areas to file
	if ($param{area} eq("yes")){
		print "The top leaflet lipid areas will be printed to ".$output_filename.".top_areas.dat\n";
		open(OUTFILE4, ">", $output_filename.".top_areas.dat");
		print OUTFILE4 @top_areas_sorted;
		close(OUTFILE4);
	}
	
	# Print bottom leaflet areas to file
	if ($param{area} eq("yes")){
		print "The bottom leaflet lipid areas will be printed to ".$output_filename.".bottom_areas.dat\n";
		open(OUTFILE5, ">", $output_filename.".bottom_areas.dat");
		print OUTFILE5 @bottom_areas_sorted;
		close(OUTFILE5);
	}


	return;
}


sub limits{
	my ($reference) = @_[0];
	my @array = @$reference;
	my $min = 1000000;
	my $max = -1000000;

	foreach (@array){
		my @temp = split(' ', $_);
		if ($temp[@_[1]-1] >= $max){$max = $temp[@_[1]-1];}
		if ($temp[@_[1]-1] <= $min){$min = $temp[@_[1]-1];}
	}

	return ($max, $min);
}


sub makegrid{
	my $grid_x = @_[0];
	my @grid_x = @$grid_x;
	my $grid_y = @_[1];
	my @grid_y = @$grid_y;
	my @grid_final;

	foreach my $n (@grid_x){
		chomp($n);
		foreach my $m (@grid_y){
			chomp($m);
			push(@grid_final, $n." ".$m."\n");
		}
	}

	return @grid_final;
}


sub convert_to_matrix{
	my ($data_name) = @_[0];
	my @data = @$data_name;
	my $grid_x = @_[1];
	my $grid_y = @_[2];
	my $prefix = @_[3];
	my $suffix = @_[4];

	open (OUT, ">", $prefix.$suffix.".dat");

	for (my $y=0; $y<$grid_y; $y++){
		for (my $x=0; $x<$grid_x; $x++){
			my $z = shift(@data);
			printf OUT "%.3f\t", $z;
		}
		print OUT "\n";
	}

	close (OUT);
	return;
}


sub convert_to_vector{
	my ($data_name) = @_[0];
	my @data = @$data_name;
	my $grid_x = @_[1];
	my $grid_y = @_[2];
	my $prefix = @_[3];
	my $suffix = @_[4];

	open (OUT, ">", $prefix.$suffix.".dat");

	for (my $y=0; $y<$grid_y; $y++){
		for (my $x=0; $x<$grid_x; $x++){
			my $z = shift(@data);
			printf OUT "%.0f\t%.0f\t%.3f\n", $x, $y, $z;
		}
	}

	close (OUT);
	return;
}



