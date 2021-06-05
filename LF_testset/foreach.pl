use File::Copy;

opendir (my $dir, "."); 
my @files = readdir($dir); 
my @pbdid; 
foreach (@files) {
	if ($_=~/^[a-z0-9]{4}$/ && -d $_) {
		push @pdbid, $_;
	}
}
close $dir;


my @names = ("xconf2-500-20.sdf", "xconf2-100-20.sdf");

foreach $pdb(@pdbid) {
	
	print "$pdb\n";
	
	# system "mkdir $pdb/actives";
	# opendir (my $dir, "actives/$pdb");
	# my @ligands = readdir($dir);
	# foreach (@ligands) {
	#	if ($_=~/\.mol/) {
	#		copy ("actives/$pdb/$_", "$pdb/actives/$_");
	#	}
	#}
	#close $dir;
	
	opendir (my $dir, "$pdb/actives"); 
	my @ligands = readdir($dir);
	open my $fo, ">$pdb/actives.sdf";
	foreach $li(@ligands) {
		if ($li =~ /\.mol$/) {
		# doing this hard way because some of ligands have trailing $$$$, others don't 
			open my $fi, "$pdb/actives/$li";
			foreach (<$fi>) {
				unless ($_=~/\$\$\$\$/) {
					print $fo $_;
				}
			}
			print $fo '$$$$'."\n";
			close $fi;
		}
	}
	close $fo;
	
	
	
	
	# system "xedmin -v $pdb/xed-ligand.sdf > $pdb/xmin-ligand.sdf";
	# system "babel  -isdf $pdb/xmin-ligand.sdf -omol $pdb/xmin-ligand.mol";
	# system "babel  -imol $pdb/xmin-ligand.mol -opdb $pdb/xmin-ligand.pdb";
	
	# system "babel -isdf $pdb/orig-ligand_opt.sdf -omol $pdb/orig-ligand_opt.mol";
	# system "babel -isdf $pdb/orig-ligand.sdf -omol $pdb/orig-ligand.mol";
	
	
	# foreach $name(@names) {
		# rename "$pdb/$name", "$pdb/orig-$name";
	# }
	
	#unless (-e "$pdb/xconf3-100.sdf")  { system "perl gen_confs.pl $pdb"; }
	#unless (-e "$pdb/xconf3-100.pdb")  { system "align -t $pdb/xmin-ligand.mol -i $pdb/xconf3-100.sdf -o $pdb/xconf3-100.pdb -olog $pdb/xconf3-100-fit.log"; }
	# unless (-e "$pdb/xconf3-100-fit.pdb")  { system "align -t $pdb/orig-ligand.mol -i $pdb/xconf3-100.sdf -o $pdb/xconf3-100-fit.pdb -olog $pdb/xconf3-100-fit2.log -sort"; }
	
	#system "mkdir update\\$pdb";
	#foreach $name(@names) {
	#	copy ("$pdb/$name", "update/$pdb/$name");
	#}
	
	#---- Metals
	# system "../../../bin/metal_1912a10.exe -og test.bin -grid -mm $pdb/bm-protein.pdb -lr $pdb/bm-reference.pdb -metal reference only";
	# open fi, "metals.log";
	# $tmp = <fi>;
	# open fo, ">>all-metals.log";
	# foreach (<fi>) { print fo "$pdb\t$_"; }
	# close fi;
	# close fo;
	# unlink "metals.log";
	
	#---- Grids
	# foreach $name(@names) {
		# unlink "$pdb/$name";
	# }
	# unlink "$pdb/bm-grid-1.bin";
	
	#-------------------
	# copy ("$pdb/xmin-ligand.mol", "$pdb/xclean-ligand.mol");
	# copy ("$pdb/bm-protein.pdb",  "$pdb/clean-protein.pdb");
	
	
	#---- Preparing stripped model
	# system "../../../bin/delhydrogens.exe $pdb/bm-protein.pdb $pdb/bm-protein-noh.pdb";
	# system "build_model -f $pdb/bm-protein-noh.pdb -omm $pdb/bm-protein-noreference.pdb -noligand";
	
	#---- Analyzing dihedrals
	#system "../../../bin/dih_stats_03.exe $pdb/xmin-ligand.mol xmin-dihedrals-03.tsv";
	#system "../../../bin/dih_stats_04.exe $pdb/xmin-ligand.mol xmin-dihedrals-04.tsv";
	#system "../../../bin/dih_stats_05.exe $pdb/xmin-ligand.mol xmin-dihedrals-05.tsv";
	
	#---- Analyzing pi-stacking
	#system "../../../bin/get_pistacking.exe -mm $pdb/bm-protein.pdb -li $pdb/bm-ligand.mol -o pistacking.tsv -a";
	
}

