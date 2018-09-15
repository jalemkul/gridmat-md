#!/usr/bin/perl

#
# convert_to_gnuplot.pl - convert the single-column matrix data output from GridMAT-MD
# into a format that can be read by Gnuplot
#
# For use with GridMAT-MD, version 1.0.2 and later
#
# GridMAT-MD was developed by W. J. Allen and J. A. Lemkul in the lab of D. R. Bevan
# If you use GridMAT-MD, please cite the paper quoted in the output of that script
#
# GridMAT-MD is now hosted at https://github.com/jalemkul/gridmat-md 
#

use strict;

# open the input from the command line
unless (scalar(@ARGV) == 3)
{
	die "Usage: $0 <input.dat> <x-grid points> <y-grid points>\n";
}

open(IN, $ARGV[0]) || die "Cannot open input: $!\n";
chomp(my @input = <IN>);
close(IN);

my $name = $ARGV[0];
my @name_array = split('\.', $name);
my $short_name = $name_array[0];

# get grid point information
my $x_grid = $ARGV[1];
my $y_grid = $ARGV[2];

# Data format: x y Z-value
# i.e., 1 1 3.295; 2 1 3.320; etc.
# nested counters necessary

# strip out unnecessary header line(s) from @input
my @input_clean;

foreach $_ (@input)
{
    unless ($_ =~ /^#/)
    {
        push(@input_clean, $_);
    }
}

# gnuplot matrix format:
# z11 z12 z13 ...
# z21 z22 z23 ...
# ...
# where each point is a Z-value for an (x,y) pair

my $out_name = $short_name."_gnuplot.dat";

print "\nOutput file: $out_name\n\n";

open(OUT, ">$out_name") || die "Cannot open output: $!\n";
print OUT "# Matrix data: Z-values for each (x,y) pair\n";

for (my $y = 1; $y < $y_grid; $y++)
{
	for (my $x = 1; $x < $x_grid; $x++)
	{
        # remove each line sequentially = don't worry about shifting indices?
        my $tmp = shift(@input_clean);
        my @tmp2 = split(" ", $tmp);
        my $z = $tmp2[2];
        printf OUT "%.3f\t", $z;
	}
	# newline at end of each x-row
	printf OUT "\n";
}

close(OUT);

exit;
