#!/usr/bin/env perl 

use strict;
use warnings;

use LWP::Simple;
use Text::CSV_XS;


#tables
my %const = ( 
	'ENA-CHECKLIST'	=> 'ERC000011',
	'genotype'	=> 'Wild type', # we only add wild-types this time
	'BioSampleModel'	=> 'Model organism or animal',
);

my %fixed = (
	"study" => "PRJNA484134",
);

my $expdate = '2013-08-25';
my $len = 51;

my $md5 = 'mouse.md5';
my %md5sum;
my $sample_tsv = 'filereport_read_run_PRJNA484134.tsv';
my %samplemap;

# my $datefile = 'mouse.txt';
# my %expdate;

my $sheet_name = "Sheet1";
#sheet_id = "1B1j5HwMJzHmIAISvxupZA5eNQz8XecW4StLRuTkMkFk"
#url = f”https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}"
my $shared_url = "https://docs.google.com/spreadsheets/d/1B1j5HwMJzHmIAISvxupZA5eNQz8XecW4StLRuTkMkFk/edit#gid=0";
my $url = $shared_url =~ s{/edit#gid=0$}{/gviz/tq?tqx=out:csv&sheet=$sheet_name}r;
my $mousetable = 'mousetable.csv';


#input MD5
{
	open my $fh, '<', $md5
		or die "cannot open ${md5}: $!";
	while(<$fh>) {
		my ($sum, $file) = split;
		$md5sum{$file} = $sum;
	}
	close $fh;
}

#input table mapping
{
	my $csv = Text::CSV_XS->new ( { binary => 1, sep_char => "\t" } )  # should set binary attribute.
		or die "Cannot use CSV: ".Text::CSV_XS->error_diag ();
 
	open my $fh, "<:encoding(utf8)", $sample_tsv 
		or die "$sample_tsv: $!";

	my @cols = @{$csv->getline ($fh)};
	my $row = {};
	$csv->bind_columns (\@{$row}{@cols});
	while ($csv->getline ($fh)) {
		$samplemap{$row->{'library_name'}} = $row->{'sample_accession'};
	}
	close $fh;
}

#input sample table
#print get($url);
mirror($url, $mousetable);

#date
## 460M    ./357_N1/20130825064613272-60416835/BSSE_QGF_13777_130821_SN792_0264_BD2C5AACXX_3_ATCACG_L003_R1_001.fastq.gz


my $csv = Text::CSV_XS->new ( { binary => 1 } )  # should set binary attribute.
	or die "Cannot use CSV: ".Text::CSV_XS->error_diag ();

open my $fh, "<:encoding(utf8)", $mousetable 
	or die "$mousetable: $!";

my @cols = @{$csv->getline ($fh)};
my $row = {};
$csv->bind_columns (\@{$row}{@cols});

#output XMLs
open my $sx, '>', 'samples.xml'
	or die 'cannot write samples.xml';
print $sx <<'SAMHEADER';
<?xml version="1.0" encoding="UTF-8"?>
<SAMPLE_SET>
SAMHEADER
my %seensample; # keep track of samples

open my $ex, '>', 'exp.xml'
	or die 'cannot write exp.xml';
print $ex <<'EXPHEADER';
<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_SET>
EXPHEADER

open my $ux, '>', 'run.xml'
	or die 'cannot write run.xml';
print $ux <<'RUNHEADER';
<?xml version="1.0" encoding="UTF-8"?>
<RUN_SET>
RUNHEADER


while ($csv->getline ($fh)) {
	my $id = $row->{'Sample Alias'};
	my $alias = $id =~ s{_}{}r;
	my $file = $row->{'FILENAME '};

	# file
	unless(defined $md5sum{$file}) {
		print STDERR "$file not in md5sum list!\n";
		exit 1;
	}

	# batch type
	my $platform = $row->{'instrument_platform'};
	my $instr = $row->{'instrument_model'};
	my $layout = $row->{'library_layout'};	# only for PAIRED NOMINAL_LENGTH="%u"
	my $strategy = $row->{'LIBRARY_STRATEGY'};
	my $source = $row->{'LIBRARY_SOURCE'};
	my $center = $row->{'Center Name'};

	my $samid = my $expid = $id;
	my $title = my $samtitle = $row->{'Sample Alias'};
	my $samdesc = '';
	
	if (defined $samplemap{$alias}) {
		print "$file\t$md5sum{$file}\t$alias\t$samplemap{$alias}\n";
		$samid=$samplemap{$alias};
		$samdesc = "accession=\"$samid\"";
	} else {
		$samdesc = "refname=\"$samid\"";
		unless ($id =~ m{^C20_}) {
			print "missing $id\n";
			exit 1;
		}

		print "$file\t$md5sum{$file}\t$alias\tcontrol\n";

		my %attr = ( %const );
		$attr{'isolate'} = $row->{'Isolate'};
		$attr{'age'} = $row->{'Age'};
		$attr{'sex'} = $row->{'Sex'};
		$attr{'tissue'} = $row->{'Tissue'};
		$attr{'disease'} = $row->{'Disease'};
		#$attr{'genotype'} = 'Wild type';
		$attr{'sample_type'} = $row->{'Sample Alias'};

	
		print $sx <<SAMPLE_IN;
 <SAMPLE alias="$samid" center_name="$center">
  <TITLE>$samtitle</TITLE>
  <SAMPLE_NAME>
   <TAXON_ID>10090</TAXON_ID>
   <SCIENTIFIC_NAME>Mus musculus</SCIENTIFIC_NAME>
   <COMMON_NAME>house mouse</COMMON_NAME>
  </SAMPLE_NAME>
  <SAMPLE_ATTRIBUTES>
SAMPLE_IN
		while (my ($tag, $value) = each (%attr)) {
			# NOTE units mandatory, _even_ when not a number
			my $unit = ($tag =~ m{longitude|latitude}) ? "\n    <UNITS>DD</UNITS>" : '';
			print $sx <<ATTR if (defined $value);
   <SAMPLE_ATTRIBUTE>
    <TAG>$tag</TAG>
    <VALUE>$value</VALUE>$unit
   </SAMPLE_ATTRIBUTE>
ATTR
		}
		print $sx <<SAMPLE_OUT;
  </SAMPLE_ATTRIBUTES>
 </SAMPLE>
SAMPLE_OUT

	}

	# append experiment
	print $ex <<EXP;
 <EXPERIMENT alias="exp-$expid">
  <TITLE>$title</TITLE>
  <STUDY_REF accession="$fixed{'study'}"/>
  <DESIGN>
   <DESIGN_DESCRIPTION/>
   <SAMPLE_DESCRIPTOR $samdesc/>
   <LIBRARY_DESCRIPTOR>
    <LIBRARY_NAME/>
    <LIBRARY_STRATEGY>$strategy</LIBRARY_STRATEGY>
    <LIBRARY_SOURCE>$source</LIBRARY_SOURCE>
    <LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>
    <LIBRARY_LAYOUT>
     <$layout/>
    </LIBRARY_LAYOUT>
    <LIBRARY_CONSTRUCTION_PROTOCOL/>
   </LIBRARY_DESCRIPTOR>
  </DESIGN>
  <PLATFORM>
   <$platform>
    <INSTRUMENT_MODEL>$instr</INSTRUMENT_MODEL>
   </$platform>
  </PLATFORM>
 </EXPERIMENT>
EXP

	# append experiment
	print $ux <<RUN;
 <RUN alias="run-$expid" center_name="$center">
  <EXPERIMENT_REF refname="exp-$expid"/>
  <DATA_BLOCK>
   <FILES>
    <FILE filename="$file" filetype="fastq"
          checksum_method="MD5" checksum="$md5sum{$file}"/>
   </FILES>
  </DATA_BLOCK>
  <RUN_ATTRIBUTES>
   <RUN_ATTRIBUTE>
    <TAG>run date</TAG>
    <VALUE>$expdate</VALUE>
   </RUN_ATTRIBUTE>
  </RUN_ATTRIBUTES>
 </RUN>
RUN
}

# finish and close everything
close($fh);

print $sx '</SAMPLE_SET>';
close($sx);

print $ex '</EXPERIMENT_SET>';
close($ex);

print $ux '</RUN_SET>';
close($ux);

exit 0;
