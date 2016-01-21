#!/usr/bin/perl
use warnings;
use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use POSIX ":sys_wait_h";

#define global variables
our $CARETAKER = "caim\@usc.edu";
our $SERVER_DIRECTORY = "/var/www/html/w4cseq";
our $HTML_DIRECTORY = "${SERVER_DIRECTORY}/html";
our $BIN_DIRECTORY = "${SERVER_DIRECTORY}/bin";
our $LIB_DIRECTORY = "${SERVER_DIRECTORY}/lib";
our $WORK_DIRECTORY = "${SERVER_DIRECTORY}/work";
our $SUBMISSION_ID_FILENAME = "${WORK_DIRECTORY}/submission_id"; #file stores the submission id

##########################################################################

$SIG{CHLD} = 'IGNORE';

my $q = new CGI;

#THE FOLLOWING LINE IS REQUIRED FOR ANY OTHER OUTPUT SUCH AS ERROR OUTPUT FOR FEEDBACK OUTPUT
print $q->header ("text/html");
my $cgi_error = $q->cgi_error;
if ($cgi_error) {
	if ($cgi_error eq "413 Request entity too large") {
		display_error ("Unable to process the information you supplied: $cgi_error! We can not handle files over 200Mb (gz/zip is okay)", $CARETAKER);
	} else {
		display_error ("Unable to process the information you supplied: $cgi_error!", $CARETAKER);
	}
}

#parse the web form at <http://w4cseq.usc.edu>
my $ip = $q->param('ip');
my $host = $q->param('host');
my $email = $q->param('email');
if ($email and $email !~ /^(\w|\-|\_|\.)+\@((\w|\-|\_)+\.)+[a-zA-Z]{2,}$/) {
	display_error("Please enter a valid email address", $CARETAKER);
}

#my $frag_method = $q->param('frag_method');
my $ref = $q->param('ref');
my $query1 = $q->param('query1');
my $query1_fh = $q->upload ('query1');
my $query2 = $q->param('query2');
my $query2_fh = $q->upload ('query2');
my $unzip="NA";
my $chipdata = "NA";
my $submission_time = scalar (localtime);
my $submission_id;
my $warning_message = '';
my $weblink;


my $chipinput = $q->param('upload_file0');
my $chipinput_fh = $q->upload('upload_file0');
my $chip1 = $q->param('upload_file1');
my $chip1_fh = $q->upload('upload_file1');
my $chip2 = $q->param('upload_file2');
my $chip2_fh = $q->upload('upload_file2');
my $chip3 = $q->param('upload_file3');
my $chip3_fh = $q->upload('upload_file3');
my $chip4 = $q->param('upload_file4');
my $chip4_fh = $q->upload('upload_file4');
my $chip5 = $q->param('upload_file5');
my $chip5_fh = $q->upload('upload_file5');
my $chip6 = $q->param('upload_file6');
my $chip6_fh = $q->upload('upload_file6');
my $chip7 = $q->param('upload_file7');
my $chip7_fh = $q->upload('upload_file7');
my $chip8 = $q->param('upload_file8');
my $chip8_fh = $q->upload('upload_file8');
my $chip9 = $q->param('upload_file9');
my $chip9_fh = $q->upload('upload_file9');
my $chip10 = $q->param('upload_file10');
my $chip10_fh = $q->upload('upload_file10');

my $chip1_name = $q->param('upload_name1');
my $chip2_name = $q->param('upload_name2');
my $chip3_name = $q->param('upload_name3');
my $chip4_name = $q->param('upload_name4');
my $chip5_name = $q->param('upload_name5');
my $chip6_name = $q->param('upload_name6');
my $chip7_name = $q->param('upload_name7');
my $chip8_name = $q->param('upload_name8');
my $chip9_name = $q->param('upload_name9');
my $chip10_name = $q->param('upload_name10');



prepareWorkDirectory ();

generateFeedback ($submission_id, $submission_time, $weblink, $warning_message);

executeProgram ();

sub executeProgram {

        if ($email eq 'caim@usc.edu') {
        	1;
        } else {
		my $pid = fork ();
		if ($pid) {	#parent process where pid=child
			1;
		} elsif ($pid == 0) {	#child process
			open STDIN,  '<', '/dev/null';		#these three lines are absolutely necessary!!!
			open STDOUT, '>', '/dev/null';		# so that the parent do not wait for child to close input/output (this is a Perl CGI issue)
			open STDERR, '>&STDOUT';
			exec("$BIN_DIRECTORY/control_w4cseq_sonication.pl $submission_id 2> $WORK_DIRECTORY/$submission_id/tempfile"); #excute control file
			exit;
		} else {
			confess "Fork process failed:$!\n";
		}
	}
}
	
sub prepareWorkDirectory {
        -d $WORK_DIRECTORY or confess "Error: work directory $WORK_DIRECTORY does not exist";
        -f $SUBMISSION_ID_FILENAME or confess "Error: submission_id file does not exist in work directory $WORK_DIRECTORY";

        open (SUBMISSION_ID, $SUBMISSION_ID_FILENAME) or confess "Error: cannot open submission_id file $SUBMISSION_ID_FILENAME in work directory $WORK_DIRECTORY: $!";
        flock SUBMISSION_ID, 1;
        $submission_id = <SUBMISSION_ID>;
        flock SUBMISSION_ID, 8;
        close (SUBMISSION_ID);
        $submission_id++;
        open (SUBMISSION_ID, ">$SUBMISSION_ID_FILENAME") or confess "Error: cannot write submission_id file in work directory $WORK_DIRECTORY: $!";
        flock SUBMISSION_ID, 2;
        print SUBMISSION_ID $submission_id;
        flock SUBMISSION_ID, 8;
        close (SUBMISSION_ID);
        
	#Xiao added the following line to generate the result link ahead of running annovar.
	my $maxLenth=16;
        my @a = (0..9,'a'..'z','A'..'Z','-','_');
        my $password = join '', map { $a[int rand @a] } 0..($maxLenth-1);        
        $weblink = qq (http://w4cseq.usc.edu/done/$submission_id/$password/index.html) ;
               
	mkdir ("$WORK_DIRECTORY/$submission_id") or confess "Error: cannot generate submission directory for submission id $submission_id: $!";
	chmod 0777, "$WORK_DIRECTORY/$submission_id" or confess "Error: unable to set the permission of directories: $!";
        
	my $orig_file1 = "$WORK_DIRECTORY/$submission_id/query1.fq";
        open (FASTQ, ">$orig_file1") or confess "Error: cannot write query1.fq file: $!";
        while (<$query1_fh>) {
                print FASTQ;
        }
        close (FASTQ);

        my $orig_file2 = "$WORK_DIRECTORY/$submission_id/query2.fq";
        open (FASTQ, ">$orig_file2") or confess "Error: cannot write query2.fq file: $!";
        while (<$query2_fh>) {
                print FASTQ;
        }
        close (FASTQ);

#generate ChIP_name file
        my $chip_names = "$WORK_DIRECTORY/$submission_id/chip_name.txt";
        open (CHIP_NAME, ">$chip_names") or confess "Error: cannot write chip_name.txt file: $!";

	if($chipinput ne "") {
		my $orig_file3 = "$WORK_DIRECTORY/$submission_id/chipc.bed";

		open (FASTQ, ">$orig_file3") or confess "Error: cannot write chipc.bed file: $!";
		while (<$chipinput_fh>) {
			print FASTQ;
		}
		close (FASTQ);
	}

        if($chip1 ne "") {
		my $orig_file4 = "$WORK_DIRECTORY/$submission_id/chip1.bed";
		print CHIP_NAME $chip1_name,"\n";
		open (FASTQ, ">$orig_file4") or confess "Error: cannot write chip1.bed file: $!";
		while (<$chip1_fh>) {
			print FASTQ;
		}
		close (FASTQ);
        }

        if($chip2 ne "") {
		my $orig_file5 = "$WORK_DIRECTORY/$submission_id/chip2.bed";
		print CHIP_NAME $chip2_name,"\n";
		open (FASTQ, ">$orig_file5") or confess "Error: cannot write chip2.bed file: $!";
		while (<$chip2_fh>) {
			print FASTQ;
		}
		close (FASTQ);
        }

        if($chip3 ne "") {
		my $orig_file6 = "$WORK_DIRECTORY/$submission_id/chip3.bed";
		print CHIP_NAME $chip3_name,"\n";
		open (FASTQ, ">$orig_file6") or confess "Error: cannot write chip3.bed file: $!";
		while (<$chip3_fh>) {
			print FASTQ;
		}
		close (FASTQ);
        }

        if($chip4 ne "") {
		my $orig_file7 = "$WORK_DIRECTORY/$submission_id/chip4.bed";
		print CHIP_NAME $chip4_name,"\n";
		open (FASTQ, ">$orig_file7") or confess "Error: cannot write chip4.bed file: $!";
		while (<$chip4_fh>) {
			print FASTQ;
		}
		close (FASTQ);
        }

        if($chip5 ne "") {
		my $orig_file8 = "$WORK_DIRECTORY/$submission_id/chip5.bed";
		print CHIP_NAME $chip5_name,"\n";
		open (FASTQ, ">$orig_file8") or confess "Error: cannot write chip5.bed file: $!";
		while (<$chip5_fh>) {
			print FASTQ;
		}
		close (FASTQ);
        }

        if($chip6 ne "") {
		my $orig_file9 = "$WORK_DIRECTORY/$submission_id/chip6.bed";
		print CHIP_NAME $chip6_name,"\n";
		open (FASTQ, ">$orig_file9") or confess "Error: cannot write chip6.bed file: $!";
		while (<$chip6_fh>) {
			print FASTQ;
		}
		close (FASTQ);
        }

        if($chip7 ne "") {
		my $orig_file10 = "$WORK_DIRECTORY/$submission_id/chip7.bed";
		print CHIP_NAME $chip7_name,"\n";
		open (FASTQ, ">$orig_file10") or confess "Error: cannot write chip7.bed file: $!";
		while (<$chip7_fh>) {
			print FASTQ;
		}
		close (FASTQ);
        }

        if($chip8 ne "") {
		my $orig_file11 = "$WORK_DIRECTORY/$submission_id/chip8.bed";
		print CHIP_NAME $chip8_name,"\n";
		open (FASTQ, ">$orig_file11") or confess "Error: cannot write chip8.bed file: $!";
		while (<$chip8_fh>) {
			print FASTQ;
		}
		close (FASTQ);
        }

        if($chip9 ne "") {
		my $orig_file12 = "$WORK_DIRECTORY/$submission_id/chip9.bed";
		print CHIP_NAME $chip9_name,"\n";
		open (FASTQ, ">$orig_file12") or confess "Error: cannot write chip9.bed file: $!";
		while (<$chip9_fh>) {
			print FASTQ;
		}
		close (FASTQ);
        }

        if($chip10 ne "") {
		my $orig_file13 = "$WORK_DIRECTORY/$submission_id/chip10.bed";
		print CHIP_NAME $chip10_name,"\n";
		open (FASTQ, ">$orig_file13") or confess "Error: cannot write chip10.bed file: $!";
		while (<$chip10_fh>) {
			print FASTQ;
		}
		close (FASTQ);
        }
	close (CHIP_NAME);

	
	if (-f $chip_names and -z _) {
		$chipdata = "no";
	}
	else {
		$chipdata= "yes";
	}	

        $unzip="yes" if $query1 !~ /\.gz$/;
        $unzip="no" if $query1 =~ /\.gz$/;

  	open (INFO, ">$WORK_DIRECTORY/$submission_id/info") or confess "Error: cannot write info file: $!";


        print INFO "email=$email\nsubmission_time=$submission_time\npassword=$password\nref=$ref\nunzip=$unzip\nquery1=$query1\nquery2=$query2";
        my $extend = $q->param('extend');
        my $bait_chr = $q->param('baitchr');
        my $bait_start = $q->param('baitstart');
        my $bait_end = $q->param('baitend');
        #my $size = $q->param('size');
	# Fan only sets up bin size for trans chromosome, so cis_bin and cis_window is not captured at all
	my $size = $q->param('trans_bin');
	
	# following three lines added by Mingyang
	my $size_inter = $q->param('trans_bin');
        my $size_intra = $q->param('cis_bin');
        my $window_intra = $q->param('cis_window');
	# end
	
        print INFO "\nextend=$extend\nbait_chr=$bait_chr\nbait_start=$bait_start\nbait_end=$bait_end\nsize=$size\nsize_inter=$size_inter\nsize_intra=$size_intra\nwindow_intra=$window_intra\nchipdata=$chipdata\n";


        close (INFO);


        if (-e "$HTML_DIRECTORY/done/$submission_id") {
        	system ("rm -rf $HTML_DIRECTORY/done/$submission_id") and confess "can not remove submission $submission_id folder";
        }
        mkdir ("$HTML_DIRECTORY/done/$submission_id") or confess "can not create submission $submission_id folder"; 
        mkdir ("$HTML_DIRECTORY/done/$submission_id/$password") or confess "can not remove submission $submission_id password $password folder";
	
        open (WAIT, ">$HTML_DIRECTORY/done/$submission_id/$password/index.html") or confess "can not create file index.html";
        print WAIT "<html><META HTTP-EQUIV=refresh CONTENT=60><p>Your submission is being processed and will be available at this page after computation is done. This page will refresh every 60 seconds. </p></html>";
        close (WAIT);
        system("chmod 777 $HTML_DIRECTORY/done/$submission_id");
        system("chmod 777 $HTML_DIRECTORY/done/$submission_id/$password");
        system("chmod 777 $HTML_DIRECTORY/done/$submission_id/$password/index.html");
}


sub generateFeedback {
	my ($submission_id, $submission_time, $weblink, $warning_message) = @_;

	
	my $submission_summary = <<SUMMARY;
<h1> Submission received </h1>
<hr>
<p>Your submission ID <b>$submission_id</b> has been received by us at <b>$submission_time</b>. </p><p>The results will be generated at <a href="$weblink"><b>$weblink</b></a> after the computation is done.</p>
<P>$warning_message</p>
SUMMARY

	print $submission_summary;
}

sub display_error
{
	my ($error_message, $email_address) = @_;
	
	my $submission_summary = "<h3> ERROR: $error_message </h3><hr><p> If this is not the expected result, please notify <var><a href=\"mailto:$email_address\">$email_address</a></var> of this error. Thank you.</p>\n";
	print $submission_summary;
	exit(0);
}


