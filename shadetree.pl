#!/usr/bin/perl

use warnings;
use strict;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::SearchIO;
use Data::Dumper;
use Bio::DB::GenBank;
use Bio::SeqIO;
use Bio::Annotation::Collection;
use Cwd;
use Getopt::Long;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::SeqIO;

## PROGRAMMING FOR BIOLOGY 2014 @ CSHL ##
## shadeTree team!! ##

  # Andrii Gryganskyi
  # Elisabeth Hehenberger
  # Iker Irisarri
  # Ryan Kepler
  # Patrick McLaughlin


#### READ CONFIG FILE ####
# IN: configuration file
# OUT: %muscle_configs hash, which contains all the info from the configuration file (not only for muscle!!)
# OUT: sends settings to blast (@blast_params) and muscle (@params)

my $config = shift @ARGV;
my %muscle_configs = ();
open(IN, '<', $config) or die "cannot open $config $!\n";
my ($ref_aln);
while (my $line = <IN>) {
    chomp $line;
    next if $line =~ /^#/;
	next if $line =~ /^\s*$/;
    my ($tag, $value) = split (/:\s+/, $line);
    $muscle_configs{$tag} =  $value;
}


if (!exists $muscle_configs{ref_aln} or $muscle_configs{ref_aln} eq ''){
    die "Please provide a reference aln\n";
}
if (!-e $muscle_configs{ref_aln} ){
    die "File does not exisit  $muscle_configs{ref_aln}\n";
}
#if (!exists $muscle_configs{params} or $muscle_configs{params} eq ''){                
#  die "Please provide a parameters for muscle run\n";                                 
#}                                                                                     
if (!exists $muscle_configs{infile} or $muscle_configs{infile} eq ''){
    die "Please provide a infile name\n";
}
if (!exists $muscle_configs{outfile} or $muscle_configs{outfile} eq ''){
    die "Please provide a outfile name\n";
}


#### GET PARAMS FOR BLAST ####
my @blast_params;
my @blast_param_keys = grep (/blast_/ , keys %muscle_configs);
foreach my $blast_param (@blast_param_keys){
    next if $blast_param eq 'blast_outfile';
    my $blast_param_val = $muscle_configs{$blast_param};
    $blast_param =~ s/blast_//;
    # get second part of blast_params: db and/or remote                                                             \
    push (@blast_params, "-$blast_param" => $blast_param_val );
}


#print join ("\n",@blast_params),"\n";


#### RUN REMOTE BLAST  ####
# IN: blast seed (fasta formatted file)
# OUT: array of accessions from significant hits @NCBI

my $cut_evalue = $muscle_configs{evalue};
print "\n\te-value cutoff: $cut_evalue\n";

print STDOUT "\nStarting BLAST search...\n";

#my $db = $muscle_configs{blast_db_name};

# create blast factory 
#run blast against ncbi server with remote = 1 or against local database by giving db name and path in -db_name
my $blast_factory = Bio::Tools::Run::StandAloneBlastPlus->new(@blast_params);
       
# run blast and store results in object $blast_result
my $blast_run = $blast_factory->blastn( -query => $muscle_configs{'query_file'},
					-outfile => $muscle_configs{'blast_outfile'},
					-method_args => [ '-num_alignments' => 20,
							  '-num_descriptions' => 20,
							  '-num_threads' => 8],
    );


my $in = new Bio::SearchIO(-format => 'blast', 
                           -file   => $muscle_configs{'blast_outfile'});

print STDOUT  "\nBLAST search done!\n",
              "\nParsing BLAST results...\n";

my %blast_results;
my $count = 0;
my $query_length;

while( my $result = $in->next_result ) {
    while( my $hit = $result->next_hit ) {
	
	# get query length of first sequence in the file
	unless (defined $query_length) {
	    $query_length = $result->query_length;
	}
	
	while( my $hsp = $hit->next_hsp ) {
	    if( $hsp->evalue < $cut_evalue ) {
		#my $query = $result->query_name;
		# get only the accession from full hit name
		my ($gb, $hit_acc, $etc) = split (/\|/, $hit->name);
		#my $evalue = $hsp->evalue;
                #my $identity = $hsp->percent_identity;
                #my $seq = $hsp->hit_string;
                # results are stored in 2D hash  
		#Save only new hits
		if (!exists $blast_results{$hit_acc} ) {
		    $blast_results{$hit_acc}{'evalue'} = $hsp->evalue;
		    $blast_results{$hit_acc}{'identity'} = $hsp->percent_identity;
		    $blast_results{$hit_acc}{'hit_aln'} = $hsp->hit_string;
		    $count++;
		} else {
		    #print "$hit_acc already present\n"
		}
                #print Dumper \%blast_results;
	    }
	}
    }
}

# this piece of code returns an array with accession numbers from significant blast results
my @all_accessions;
foreach my $key (keys %blast_results) {
    push (@all_accessions, $key);
}

#### GET STUFF FROM GENBANK  ####
# IN: array of accessions from significant hits @NCBI  
# OUT: fasta containing all sequences from significant hits
# OUT: tab-delimited file with metadata

print STDOUT "\nParsing BLAST results done!\n",
             "\n\t$count non-redundant significant hits found\n",
             "\nDownloading data from NCBI...\n",
             "\nSavig metadata to shadetree_meta.txt\n\n\t&\n",
             "\nFormatting fasta file from GenBank hits...\n";

#my @all_accessions = qw(FJ554191.1 GQ219901.1 AY704761.1);
my $gb = Bio::DB::GenBank->new(); # connects to genbank

my @org_name;
my @isol_source;
my @country;
my @date;
my %meta_data;
my $feature_object;
my %sequences;

open (OUT, '>', "shadetree.fasta") or die "can't open shadetree.fasta\n";
open (OUT2, '>', "shadetree_meta.txt") or die "can't open shadetree_meta\n";

print OUT2 "Sequence_Name\tAccession\tOrganism\tIsolation_Source\tCountry\tCollection_Date\tReference\n";

foreach my $accession (@all_accessions) {
	
        my $seq_object = $gb->get_Seq_by_version($accession); # create sequence object for each accession 

	my $annot_collection = $seq_object->annotation; # create annotation collection object to access the annotation object (everything that's not in SeqFeature objects) 
	my @annotations = $annot_collection->get_Annotations("reference"); # get the array "reference" from the annotation collection
#	print $annotations[0]{"authors"}, ",", $annotations[0]{"title"}, ",",  $annotations[0]{"location"}, "\n"; # retrieve first element of array (first reference) and value corresponding to the given key from the hash ref stored in the array		

	@org_name = "NA"; # reset variable to "NA"
	@isol_source = "NA";
	@country = "NA";
	@date = "NA";

        foreach $feature_object ($seq_object->get_SeqFeatures) { # get_SeqFeatures creates list (=array) of feature objects
        	my $primary_tag = $feature_object->primary_tag; # returns primary tag (=string) of features
		if ($primary_tag ne "source") { # skip all primary_tags that are not "source"
		next;
		}
		else {
			unless ($feature_object->has_tag("organism")==0) { # unless the tag of interest is not present (FALSE)
				@org_name = $feature_object->get_tag_values('organism'); # get the tag value (returned as a list)
				$meta_data{"organism"} = \@org_name; # store reference for the array as value for the key "organism" in %hash
			}
			unless ($feature_object->has_tag("isolation_source")==0) {
				@isol_source = $feature_object->get_tag_values("isolation_source");
                                $meta_data{"isolation_source"} = \@isol_source;
			}
			unless ($feature_object->has_tag("country")==0) {
                        	@country = $feature_object->get_tag_values("country");
			        $meta_data{"country"} = \@country;
			}	
			 unless ($feature_object->has_tag("collection_date")==0) {
                                @date = $feature_object->get_tag_values("collection_date");
                                $meta_data{"collection_date"} = \@date; 
                        }    		
		}
	}

	#my @acc_array = @org_name;
	map { $_ =~ s/\s+/_/g } @org_name; # map accesses every element in array -> becomes assigned to $_ variable -> substitute inside element
	#map { $_ =~ s/$/\@$accession/ } @acc_array;

	my $seq_name = join("", @org_name) . "\@" . $accession;
	
	my $seq = $seq_object->seq;
	$sequences{$accession}{"seq_name"} = $seq_name;
	
	if (length $seq < (2 * $query_length)) {
	    $sequences{$accession}{"sequence"} = $seq_object->seq;
	} else {
	    my $hit_aln = $blast_results{$accession}{'hit_aln'};
	    $sequences{$accession}{"sequence"} = $hit_aln;
	    print "\they, I am $seq_name and too long!! Replacing seq with hit string\n"; 
	}

	print OUT2 "$seq_name\t"; 	
	print OUT2 "$accession\t";
	print OUT2 join("", @org_name), "\t";
        print OUT2 join("", @isol_source), "\t";
        print OUT2 join("", @country), "\t";
        print OUT2 join("", @date), "\t";
	print OUT2 $annotations[0]{"authors"}, " ", $annotations[0]{"title"}, ". ",  $annotations[0]{"location"}, "\n";

}

foreach my $key (keys %sequences) {
	print OUT ">", $sequences{$key}{"seq_name"}, "\n";
	print OUT $sequences{$key}{"sequence"}, "\n";
}

close (OUT);
close (OUT2);

print STDOUT "\nDownloading NCBI data done!\n",
             "\nStarting Muscle alignment...\n" ;

#### START MULTIPLE SEQUENCE ALIGNMENT ####
# IN: 
# OUT:
 
#### GET PARAMS FOR MUSCLE ####                                                                                     

my @params;
my @param_keys = grep (/param_/ , keys %muscle_configs);
foreach my $param (@param_keys){
    my $param_val = $muscle_configs{$param};
    $param =~ s/param_//;
    push @params , ("-$param" => $param_val );
}

#print join ('--',@params),"\n";

#### RUN MUSCLE ####

print "\nGood to Go!!\n",
      "\nConcatenating new blast hits with reference alignment\n";

# concatenate new fasta file with reference alignment
`cat $muscle_configs{ref_aln} $muscle_configs{infile} > cat.fa`;


# The "factory" created for the alignment with the parameters above     
#my $aln_factory = Bio::Tools::Run::Alignment::Muscle->new($muscle_configs{params});

my $aln_factory = Bio::Tools::Run::Alignment::Muscle->new(@params);

# The input file with protein sequences in FASTA format
my $str = Bio::SeqIO->new(-file => "cat.fa", '-format' => 'Fasta');

# The output file where the alignment will be writen - in clustalw format     
my $out = Bio::AlignIO->new(-file => ">$muscle_configs{outfile}", -format => 'fasta');

# First creates an array with all sequences to be alignment
my @allseq_array =();

while ( my $seq = $str->next_seq() ) {
    push (@allseq_array, $seq);
    print "Added sequence ".$seq->id." to the array\n";
}

# Then align the sequences in the array using the factory created before
my $aln = $aln_factory->align(\@allseq_array);

$out->write_aln($aln);

print "\nTrimming multiple sequence alignment...\n";

my $trimal_cmd = "trimal -in $muscle_configs{outfile} -out $muscle_configs{trimalout} -fasta -gt 0.9 -cons 60";  #trimal command line    
my $call = system $trimal_cmd;

die "Trimal error $!\n" if $call;

print "\nStarting ML search...\n";

my $raxml_directory = $muscle_configs{raxml_dir};
my $raxml_mkdir = system "mkdir $muscle_configs{raxml_dir}; chmod a+w $muscle_configs{raxml_dir}";
die "Unable to make $muscle_configs{raxml_dir} $!\n" if $raxml_mkdir;

my $path_name = getcwd;
my $raxml_cmd = "raxmlHPC-PTHREADS -T $muscle_configs{raxml_cores} -p 88998 -x 99889 -f a -# $muscle_configs{raxml_boots} -m $muscle_configs{raxml_model} --no-bfgs -s $muscle_configs{trimalout} -n $muscle_configs{run_name} -w $path_name/$muscle_configs{raxml_dir} ";

print $raxml_cmd;

my $raxml_call = system $raxml_cmd;

die "your raxml command call is bad!!! $!\n" if $raxml_call;
