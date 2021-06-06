#!/usr/bin/env perl 

use strict;
use warnings;

#tables
my %wwtp = (
	10 => [ 'ZH - Zürich - ARA Werdhölzli',	47.4007,	8.4810 ],
	12 => [ 'VD - Lausanne - STEP Vidy',	46.5202,	6.5929 ],
	24 => [ 'Ski Resort', 'restricted access', 'restricted access' ], # disclosure of the ski resort is not allowed
);

my %const = ( 
	"ENA-CHECKLIST"	=> "ERC000023",
	"project name"	=> "Swiss SARS-CoV-2 Wastewater Surveillance",
	"sequencing method"	=> "NGS Illumina xxx",
	"investigation type"	=> "virus",
	"geographic location (country and/or sea)"	=> "Switzerland",
	"wastewater/sludge environmental package"	=> "wastewater/sludge",
	"environment (biome)"	=> "urban biome",
	"environment (feature)"	=> "city",
	"environment (material)"	=> "sewage",
	"sewage type"	=> "municipal wastewater treatment plant",
	"wastewater type"	=> "human waste",
);

my %fixed = (
	"study" => "PRJEB44932", # "PRJEB44932",
);

my $rx_sam = qr{(?P<samid>(?:(?P<plant>\d+)_(?P<year>20\d{2})_(?:(?:(?P<month>[01]?\d)_(?P<day>[0-3]?\d))|(?:R_(?P<repeat>\d+))))|(?:Wild_(?P<dilution>\d+)_(?P<dilrep>\d+)))};
my $rx_bat = qr{(?P<year>20\d{2})(?P<month>[01]\d)(?P<day>[0-3]?\d)_(?P<flow>\w)};


# parameter
my $listfile = shift @ARGV;
my $md5 = shift @ARGV;
my %md5sum;


#input MD5
{
	open my $fh, '<', $md5
		or die "cannot open ${md5}";
	while(<$fh>) {
		my ($sum, $file) = split;
		$md5sum{$file} = $sum;
	}
	close $fh;
}

#input TSV
open my $fh, '<', $listfile
	or die "cannot open ${listfile}";

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


while(<$fh>) {
	my ($sam, $bat, $len) = split;
	next 
		unless ($sam =~ $rx_sam);

	my $expid = $+{'samid'};
	my $title = my $samtitle = my $samid = '';
	my %attr = ( %const );

	# sample types
	if (defined $+{'dilrep'}) {
		# special case: dilution serie

		$attr{'collection date'} = '2021-03-09';
		$attr{'geographic location (latitude)'} =
		$attr{'geographic location (longitude)'} = 'not provided'; # it's an artificial mixture, it doesn't really have a sampling place
		my $dilution = $+{'dilution'};
		my $rep = $+{'dilrep'};
		$samid = $expid =~ s{_$rep$}{}r;
		$samtitle = "SARS-CoV2 Wild:B117 1:$dilution dilution wastewater sample";
		$title = "$samtitle replicate $rep";
	} elsif (defined $+{'repeat'}) {
		# special case: replicate serie

		$attr{'collection date'} = '2021-03-09';
		( $attr{'geographic location (region and locality)'},
		$attr{'geographic location (latitude)'},
		$attr{'geographic location (longitude)'}, ) = @{$wwtp{int $+{'plant'}}};
		my $rep = $+{'repeat'};
		$samid = $expid =~ s{_$rep$}{}r;
		$samtitle = "SARS-CoV2 replication serie wastewater sample";
		$title = "$samtitle replicate $rep";
	} else {
		# standard waste-water sample
		next
			unless (defined $wwtp{int $+{'plant'}});

		$attr{'collection date'} = sprintf('%04u-%02u-%02u', $+{'year'}, $+{'month'}, $+{'day'});
		( $attr{'geographic location (region and locality)'},
		$attr{'geographic location (latitude)'},
		$attr{'geographic location (longitude)'}, ) = @{$wwtp{int $+{'plant'}}};
		$samid = $expid;
		$samtitle = $title = "SARS-CoV-2 positive wastewater sample $samid";

		# HACK triplicate, within the same batch
		$expid .= '_3' if ($sam =~ m{_n});

		# append batch
		$expid .= "_$bat";
		
	}
	print "<$samid>\t$title\n";

	# batch type
	# HACK recognize the instrument based ond the flow cell id
	my $instr = 'Illumina';
	$bat =~ $rx_bat;
	my $expdate = "$+{'year'}-$+{'month'}-$+{'day'}";
	if ($+{'flow'} eq 'J') { # e.g.:_JG6YF
		$instr .= ' MiSeq';
	} elsif ($+{'flow'} eq 'H') { # e.g.: _HYW2CDRXX
		$instr .= ' NovaSeq 6000';
	}

	# file
	my $file = "${sam}_${bat}.cram";
	unless(defined $md5sum{$file}) {
		print STDERR "$file not in md5sum list!\n";
		exit 1;
	}

	# append sample.xml
	unless(defined $seensample{$samid}) {
		# NOTE only generate a single sample. Not for each repeated sequencing
		print $sx <<SAMPLE_IN;
 <SAMPLE alias="$samid" center_name="">
  <TITLE>$samtitle</TITLE>
  <SAMPLE_NAME>
   <TAXON_ID>2697049</TAXON_ID>
   <SCIENTIFIC_NAME>Severe acute respiratory syndrome coronavirus 2</SCIENTIFIC_NAME>
   <COMMON_NAME>SARS-CoV2</COMMON_NAME>
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
		++$seensample{$samid};
	}

	# append experiment
	print $ex <<EXP;
 <EXPERIMENT alias="exp-$expid">
  <TITLE>$title</TITLE>
  <STUDY_REF accession="$fixed{'study'}"/>
  <DESIGN>
   <DESIGN_DESCRIPTION/>
   <SAMPLE_DESCRIPTOR refname="$samid"/>
   <LIBRARY_DESCRIPTOR>
    <LIBRARY_NAME/>
    <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
    <LIBRARY_SOURCE>METAGENOMIC</LIBRARY_SOURCE>
    <LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>
    <LIBRARY_LAYOUT>
     <PAIRED NOMINAL_LENGTH="$len"/>
    </LIBRARY_LAYOUT>
    <LIBRARY_CONSTRUCTION_PROTOCOL/>
   </LIBRARY_DESCRIPTOR>
  </DESIGN>
  <PLATFORM>
   <ILLUMINA>
    <INSTRUMENT_MODEL>$instr</INSTRUMENT_MODEL>
    </ILLUMINA>
  </PLATFORM>
 </EXPERIMENT>
EXP

	# append experiment
	print $ux <<RUN;
 <RUN alias="run-$expid" center_name="">
  <EXPERIMENT_REF refname="exp-$expid"/>
  <DATA_BLOCK>
   <FILES>
    <FILE filename="$file" filetype="cram"
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

close($fh);

print $sx '</SAMPLE_SET>';
close($sx);

print $ex '</EXPERIMENT_SET>';
close($ex);

print $ux '</RUN_SET>';
close($ux);

exit 0;
