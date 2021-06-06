
my $nargs = @ARGV;
if ($nargs != 1) {
	die "usage: calc_charges <input>";
}
my $structure = $ARGV[0];

my $totalcharge = "";

if ($structure =~/\.pdb$/) {
	$totalcharge = totalcharge_pdb($structure);
}
elsif ($structure =~/\.mol2$/) {
	$totalcharge = totalcharge_mol2($structure);
}
print "$totalcharge";



#===============================================================

sub totalcharge_pdb {
	my $result = 0;
	open my $fi, $structure;
	my $count = 0;

	foreach (<$fi>) {
		if ($_=~/^ATOM/ || $_=~/^HETATM/) {
			$count++;
			if ($_=~/^.{77}\S+(\d+)[+]/) {
				$result += $1;
			}
			if ($_=~/^.{77}\S+(\d+)[-]/) {
				$result -= $1;
			}
		}
	}
	close $fi;
	if ($count == 0) { return ""; }
	return $result;
}

sub totalcharge_mol2 {
	
	my $result = 0;
	my $count = 0;
	open my $fi, $structure;
	foreach (<$fi>) {
		if ($_=~/\d+\s+\S+\s+[-0-9.]+\s+[-0-9.]+\s+[-0-9.]+\s+\S+\s+\d+\s+\S+\s+(\S+)/) {
			$result += $1;
			$count++;
		}
	}
	close $fi;
	if ($count == 0) { return ""; }
	return $result;
	
}