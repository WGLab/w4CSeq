#!/usr/bin/env perl
use warnings;
use strict;
use Carp;
use POSIX ":sys_wait_h"; ##zombie process

use File::Basename;

#define global variables
our $CARETAKER = "caim\@usc.edu";

our $SERVER_DIRECTORY = "/var/www/html/w4cseq/";
our $HTML_DIRECTORY = "${SERVER_DIRECTORY}/html";

our $WORK_DIRECTORY = "$SERVER_DIRECTORY/work";
our $LIB_DIRECTORY = "$SERVER_DIRECTORY/lib";
our $BIN_DIRECTORY = "$SERVER_DIRECTORY/bin";
our $SUBMISSION_ID_FILENAME = "$WORK_DIRECTORY/submission_id";
our $PID_LOCK_FILENAME = "$WORK_DIRECTORY/pid.lock";
our $zombies = 0;

#prepare PATH environmental variable
$ENV{PATH} = "$BIN_DIRECTORY:$ENV{PATH}";

#BEGINNING THE CHECKING OF LOG, ERROR AND SUBMISSION_ID FILES, RETRIEVE THE PREVIOUS FAILURE IF THERE IS ANY, SET THE CORRECT OLD_ID VALUE
my ($old_submission_id, $submission_id, $keep_temp);

#setting up a signal handler for Ctrl+C, which can be used to stop this program from running
$SIG{INT} = \&deletePidFile;
#statistic number of zombie process
$SIG{CHLD} = sub { $zombies++ };  


if (defined $ARGV[0]) {
	$keep_temp = $ARGV[1]||0;
	open (STDERR, ">$WORK_DIRECTORY/w4cseq_error_log_enzyme") or confess "Error: unable to redirect STDERR to error log file $WORK_DIRECTORY/w4cseq_error_log_enzyme";

	if ($ARGV[0] =~ m/^(\d+)$/) {
		chdir ("$WORK_DIRECTORY/$1") or confess "Error: cannot change to the desired directory $WORK_DIRECTORY/$1: $!";
		print STDERR "\n\n\n\n\nNOTICE: time=", scalar (localtime), " command=$0\n";
		processSubmission ($ARGV[0]);
	} elsif ($ARGV[0] =~ m/^(\d+)\-(\d+)$/) {
		for ($1 .. $2) {
			chdir ("$WORK_DIRECTORY/$_") or confess "Error: cannot change to the desired directory $WORK_DIRECTORY/$_: $!";
			print STDERR "\n\n\n\n\nNOTICE: time=", scalar (localtime), " command=$0\n";
			processSubmission ($_);
		}
	} elsif ($ARGV[0] eq '-h') {
		die "Usage: $0 <submission-id> [keep-temp]\nFunction: process 4CSEQ submission. When no argument is given this program will run continuously until pressing Ctrl+C\n\n" .
			"Process submission id 10, 11, .. 20\n";
	} else {
		confess "Error: <submission-id> must be a positive integer or a range of positive integers!\n\nUsage: $0 <submission-id> [keep-temp]\nFunction: process submissions\n";
	}
} else {

	#now we want to record the STDERR so that we know what the error messages are
	open (STDERR, " | tee -a loh_error_log 1>&2") or confess "Error: unable to redirect STDERR to error log file";
	print STDERR "\n\n\n\n\nNOTICE: time=", scalar (localtime), " command=$0\n";

	checkAndWritePidFile ();
	$old_submission_id = readSubmissionID ();
	my $last_processed_submission_id = checkUnprocessedSubmission ();
	if ($last_processed_submission_id <= $old_submission_id) {
		print STDERR "*"x80, "\nWARNING: submissions from $last_processed_submission_id to $old_submission_id have NOT been processed!\n";
		print STDERR "Pressing Ctrl+C then running 'control_w4cseq.pl $last_processed_submission_id-$old_submission_id' will solve this probem!\n", "*"x80, "\n";
	}
	print STDERR "old_submission_id = $old_submission_id ...\n";
	while (1) {
		$submission_id = readSubmissionID ();
		if ($submission_id > $old_submission_id) {
			for (($old_submission_id+1) .. $submission_id) {
				while (`ps aux | fgrep control_w4cseq.pl | wc -l` > 6) {
					print "Exceed the maximum number of processing submissions";
					sleep 30;
				}
				defined(my $pid1=fork()) or die "Fork process failed:$!\n";
				if ($pid1) {
					&REAPER if $zombies;
				} else {
					print STDERR "processing submission $_ ...\n";
					processSubmission ($_);  
					exit;   
				}
			}
		}
		$old_submission_id = $submission_id;
		print STDERR "wait 100 seconds before checking again ...\n";
		select (undef, undef, undef, 100);
	}
}

#this subroutine is used to delete the temporary PID file, so that we know the program exits successfully.
#if a PID file exists but the program is not running, then we know there must be something wrong and the program exits unexpectedly.
sub deletePidFile {
	-f $PID_LOCK_FILENAME or confess "Error: cannot find pid file $PID_LOCK_FILENAME";
	unlink ($PID_LOCK_FILENAME) or confess "Error: cannot delete pid file $PID_LOCK_FILENAME: $!";
	print STDERR "Successfully deleted pid file $PID_LOCK_FILENAME\nProgram exited!\n";
	exit (0);
}

#this subroutine check the existence of PID file, and generate a new one if it does not exist.
#if a PID file already exists, then it is possible that (1) another instance of the program is running (2) the previous program exits unexpectedly
sub checkAndWritePidFile {
	-f $PID_LOCK_FILENAME and confess "Error: another process may be running. check the $PID_LOCK_FILENAME file for details";
	open (PID, ">$PID_LOCK_FILENAME") or confess "Error: cannot write to pid file $PID_LOCK_FILENAME: $!";
	print PID "time=", scalar (localtime), "\n";
	print PID "host=", `hostname`;
	print PID "pid=$$\n";
	close (PID);
}

sub readSubmissionID {
	open (SUBMISSION_ID, $SUBMISSION_ID_FILENAME) or confess "Error: cannot open submission_id file $SUBMISSION_ID_FILENAME: $!";
	flock SUBMISSION_ID, 1;
	my $submission_id = <SUBMISSION_ID>;
	flock SUBMISSION_ID, 8;
	close (SUBMISSION_ID);
	return $submission_id;
}

sub checkUnprocessedSubmission {
	my $submission_id = $old_submission_id;
	while ($submission_id) {
		-f "$WORK_DIRECTORY/done/$submission_id/email" and last;#######
		$submission_id--;
	}
	$submission_id++;
	return $submission_id;
}

sub processSubmission {
        my ($id) = @_;
        my %info;
	my ($system_command);
	my ($failed_command);
	my $result_page = '';

	#enter the working directory, so that all temporary files are confined in this directory
	chdir ("$WORK_DIRECTORY/$id") or confess "Error: cannot change to the desired directory $WORK_DIRECTORY/$id: $!";

	open (INFO, "info") or warn "Error: cannot read info file $WORK_DIRECTORY/$id/info: $!" and return;
	while (<INFO>) {
		chomp;
		my ($key, $value) = split ('=', $_);
		$info{$key} = $value;
	}
	close (INFO);

	-f "email" and -s "email" and print STDERR "WARNING: skipping submission $id since it has been already executed (the $WORK_DIRECTORY/$id/email file exists)\n" and return;

        my $process_time = scalar (localtime);
        my ($email_header, $email_body, $email_tail);
	my $password = $info{password};
	
	
	$system_command = "/var/www/html/w4cseq/bin/R-3.1.2/bin/Rscript $BIN_DIRECTORY/4C_enzyme.R 1 query1.fq $info{ref} $info{target} $info{enzyme} $info{bait_chr} $info{bait_start} $info{bait_end} $info{size_inter} $info{size_intra} $info{window_intra} $id $info{unzip} $info{chipdata}> $WORK_DIRECTORY/$id/run_log.txt";


	system ($system_command) and die "cannot run system command <$system_command>\n";


        system ("cp $WORK_DIRECTORY/$id/FASTQ_FOR_MAPPING.fq /var/www/html/w4cseq/html/done/$id/$password");
        system ("cp $WORK_DIRECTORY/$id/FASTQ_FILTERED.fq /var/www/html/w4cseq/html/done/$id/$password");
        system ("cp $WORK_DIRECTORY/$id/MAPPED_BAM.bam /var/www/html/w4cseq/html/done/$id/$password");
        system ("cp $WORK_DIRECTORY/$id/MAPPED_BAM.bam.bai /var/www/html/w4cseq/html/done/$id/$password");
        system ("cp $WORK_DIRECTORY/$id/DISTAL_INTERACTION_SITES.bed /var/www/html/w4cseq/html/done/$id/$password");
	system ("cp $WORK_DIRECTORY/$id/positive_hits.bed /var/www/html/w4cseq/html/done/$id/$password");
	system ("cp $WORK_DIRECTORY/$id/SIGNIFICANT_REGIONS.bed /var/www/html/w4cseq/html/done/$id/$password");
	system ("cp $WORK_DIRECTORY/$id/UCSC_view.bed /var/www/html/w4cseq/html/done/$id/$password");

	
	system ("cp $WORK_DIRECTORY/$id/circos.pdf /var/www/html/w4cseq/html/done/$id/$password");
	system ("cp $WORK_DIRECTORY/$id/circos.png /var/www/html/w4cseq/html/done/$id/$password");
	
	system ("cp $WORK_DIRECTORY/$id/genome.pdf /var/www/html/w4cseq/html/done/$id/$password");
	system ("cp $WORK_DIRECTORY/$id/genome.png /var/www/html/w4cseq/html/done/$id/$password");
	
	system ("cp $WORK_DIRECTORY/$id/domainogram.pdf /var/www/html/w4cseq/html/done/$id/$password");
	system ("cp $WORK_DIRECTORY/$id/domainogram.png /var/www/html/w4cseq/html/done/$id/$password");
	
	system ("cp $WORK_DIRECTORY/$id/spider.pdf /var/www/html/w4cseq/html/done/$id/$password");
	system ("cp $WORK_DIRECTORY/$id/spider.png /var/www/html/w4cseq/html/done/$id/$password");
	
        system ("cp $WORK_DIRECTORY/$id/distance.pdf /var/www/html/w4cseq/html/done/$id/$password");
	system ("cp $WORK_DIRECTORY/$id/distance.png /var/www/html/w4cseq/html/done/$id/$password");
	
        system ("cp $WORK_DIRECTORY/$id/DNA_replication.pdf /var/www/html/w4cseq/html/done/$id/$password");
	system ("cp $WORK_DIRECTORY/$id/DNA_replication.png /var/www/html/w4cseq/html/done/$id/$password");
	
	system ("cp $WORK_DIRECTORY/$id/ChIP-Seq.pdf /var/www/html/w4cseq/html/done/$id/$password");
	system ("cp $WORK_DIRECTORY/$id/ChIP-Seq.png /var/www/html/w4cseq/html/done/$id/$password");
        

	#system("wc fastq_convert.fq | awk '{print \"Total number of bait-containing reads = \",$1/4}' >> summary_report.txt")
	#system("wc FASTQ_FILTERED.fq | awk '{print \"Total number of bait-containing reads with good base quality = \",$1/4}' >> summary_report.txt")
	#system("wc all_interact.bed | awk '{print \"Total number of interaction reads after removing randomly aligned reads = \",$1}' >> summary_report.txt")
	#system("wc distal_interact.bed | awk '{print \"The number of distal interaction reads after removing randomly aligned reads = \",$1}' >> summary_report.txt")
	#system("wc DISTAL_INTERACTION_SITES.bed | awk '{print \"The number of distal interacting sites = \",$1}' >> summary_report.txt")
	#system("wc SIGNIFICANT_REGIONS.bed | awk '{print \"The number of significant interacting regions = \",$1}' >> summary_report.txt")

	my $total_reads_count = `wc -l < "$WORK_DIRECTORY/$id/fastq_convert.fq"`/4;
	my $good_reads_count = `wc -l < "$WORK_DIRECTORY/$id/FASTQ_FILTERED.fq"`/4;
	my $nonrandom_interact_reads_count = `wc -l < "$WORK_DIRECTORY/$id/all_reads.bed"`;
	#my $all_interact_count = `wc -l < "$WORK_DIRECTORY/$id/all_interact.bed"`;
	#my $distal_interact_count = `wc -l < "$WORK_DIRECTORY/$id/distal_interact.bed"`;
	my $distal_sites_count = `wc -l < "$WORK_DIRECTORY/$id/DISTAL_INTERACTION_SITES.bed"`;
	my $positive_sites_count = `wc -l < "$WORK_DIRECTORY/$id/positive_hits.bed"`;
	my $signif_regions_count = `wc -l < "$WORK_DIRECTORY/$id/SIGNIFICANT_REGIONS.bed"`;
	
 	#open(RES, "$WORK_DIRECTORY/$id/summary_report.txt") or die "Error: cannot read result table file: $!\n";

	#produce the result page	
	open (HTML, ">$HTML_DIRECTORY/done/$id/$password/index.html") or die "can't write out index.html: $!";
	print HTML "
		<!DOCTYPE>
		<html lang=\"en\">
		<head>
		    <meta http-equiv=\"Content-type\" content=\"text/html; charset=utf-8\">
                    <div class=\"container\">
                        <br><br>
                    </div>
		    <title>4C RESULTS</title>
		    <link rel=\"stylesheet\" href=\"http://w4cseq.usc.edu/W4CSEQ_files/css/bootstrap.css\" type=\"text/css\" media=\"screen\" title=\"master\" charset=\"utf-8\">
		</head>
	 
                <body>
                    
                    <div class=\"container\">
                        <header class=\"row\">
			    <h2>Dear W4CSEQ user, your submission (identifier: $id) was received at $info{submission_time} and processed at $process_time.</h2>
			    <hr>
                            <h1 id=\"banner\">4C RESULTS</h1>
                        </header>
                         
                        <nav id=\"primary-navigation\" class=\"row\">
                            <div class=\"col-sm-12\">
                                <ul>
                                    <li><a href=\"#parameters\">Summary of the parameters</a></li>
                                    <li><a href=\"#metrics\">Summary of the metrics</a></li>
                                    <li><a href=\"#figures\">Figures</a></li>
                                    <li><a href=\"#files\">Files</a></li>
				    <li><a href=\"#UCSC\">View in UCSC genome browser</a></li>
                                </ul>
                            </div>
                        </nav>
                         
                        <section class=\"row\">
                            <article class=\"col-sm-12\">
                            
                                <div class=\"col-md-8\" id = \"parameters\">
                                    <h2>Summary of parameters</h2>
                                        
                                        <table class=\"table table-hover\">
                                            <thead>
                                              <tr>
                                                <th>ITEM</th>
                                                <th>VALUE</th>
                                              </tr>
                                            </thead>
                                            <tbody>
                                                <tr>
                                                    <td>FASTQ_FILE</td>
                                                    <td>$info{query1}</td>
                                                <tr>
                                            </tbody>
                                            <tbody>
                                                <tr>
                                                    <td>GENOME_BUILD</td>
                                                    <td>$info{ref}</td>
                                                <tr>
                                            </tbody>
                                            <tbody>
                                                <tr>
                                                    <td>RESTRICTION_ENZYME</td>
                                                    <td>$info{enzyme}</td>
                                                <tr>
                                            </tbody>
                                            <tbody>
                                                <tr>
                                                    <td>BAIT_SEQUENCE</td>
                                                    <td>$info{target}</td>
                                                <tr>
                                            </tbody>
                                            <tbody>
                                                <tr>
                                                    <td>BAIT_CHROMOSOME</td>
                                                    <td>$info{bait_chr}</td>
                                                <tr>
                                            </tbody>
                                            <tbody>
                                                <tr>
                                                    <td>BAIT_START</td>
                                                    <td>$info{bait_start}</td>
                                                <tr>
                                            </tbody>
                                            <tbody>
                                                <tr>
                                                    <td>BAIT_END</td>
                                                    <td>$info{bait_end}</td>
                                                <tr>
                                            </tbody>
                                            <tbody>
                                                <tr>
                                                    <td>BIN_TRANS</td>
                                                    <td>$info{size_inter}</td>
                                                <tr>
                                            </tbody>
                                            <tbody>
                                                <tr>
                                                    <td>BIN_CIS</td>
                                                    <td>$info{size_intra}</td>
                                                <tr>
                                            </tbody>
                                            <tbody>
                                                <tr>
                                                    <td>WINDOW_CIS</td>
                                                    <td>$info{window_intra}</td>
                                                <tr>
                                            </tbody>
                                        </table>
                                        
                                </div>
                               
                                <div class=\"col-md-8\" id = \"metrics\">
                                    <h2>Summary of 4C metric</h2>
                                        <table class=\"table table-hover\">
                                        <thead>
                                          <tr>
                                            <th>ITEM</th>
                                            <th>COUNT</th>
                                          </tr>
                                        </thead>
                                        <tbody>
                                            <tr>
                                                <td>Bait-containing reads</td>
                                                <td>$total_reads_count</td>
                                            <tr>
                                        </tbody>
                                        <tbody>
                                            <tr>
                                                <td>Bait-containing reads with good base quality (>=20)</td>
                                                <td>$good_reads_count</td>
                                            <tr>
                                        </tbody>
                                        <tbody>
                                            <tr>
                                                <td>Interaction reads after removing randomly aligned reads</td>
                                                <td>$nonrandom_interact_reads_count</td>
                                            <tr>
                                        </tbody>
					
                                        <tbody>
                                            <tr>
                                                <td>Interacting sites</td>
                                                <td>$distal_sites_count</td>
                                            <tr>
                                        </tbody>
					<tbody>
                                            <tr>
                                                <td>Significant interacting sites</td>
                                                <td>$positive_sites_count</td>
                                            <tr>
                                        </tbody>
                                        <tbody>
                                            <tr>
                                                <td>Significant interacting regions</td>
                                                <td>$signif_regions_count</td>
                                            <tr>
                                        </tbody>
                                    
                                    </table>
                                        
                                
                                </div>    
                                    
                                <div class=\"col-md-11\" id = \"figures\">
                                    <h2>Figures</h2>
                                    
                                    <div class=\"col-md-11\">
                                    <div class=\"col-md-8\">
                                        <img class=\"img-thumbnail\" data-src=\"holder.js/500x500/auto\" alt=\"500x500\" src=\"genome.png\"><br>
                                        <p class=\"bg-info\">Genome-wide distribution of 4C regions</p>
                                        <br><br>
                                    </div>
                                    <div class=\"col-md-3\">
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/genome.pdf\" class=\"btn btn-primary btn-large\">Download PDF</a></p>
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/genome.png\" class=\"btn btn-primary btn-large\">Download PNG</a></p>
                                    </div>
                                    </div>
                                    
				    <div class=\"col-md-11\">
                                    <div class=\"col-md-8\">
                                        <img class=\"img-thumbnail\" data-src=\"holder.js/500x500/auto\" alt=\"500x500\" src=\"domainogram.png\"><br>
                                        <p class=\"bg-info\">Domainogram of intra-chromosomal interactions</p>
                                        <br><br>
                                    </div>
                                    <div class=\"col-md-3\">
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/domainogram.pdf\" class=\"btn btn-primary btn-large\">Download PDF</a></p>
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/domainogram.png\" class=\"btn btn-primary btn-large\">Download PNG</a></p>
                                    </div>
                                    </div>
				    
				    <div class=\"col-md-11\">
                                    <div class=\"col-md-8\">
                                        <img class=\"img-thumbnail\" data-src=\"holder.js/500x500/auto\" alt=\"500x500\" src=\"spider.png\"><br>
                                        <p class=\"bg-info\">Spider plot of intra-chromosomal interactions</p>
                                        <br><br>
                                    </div>
                                    <div class=\"col-md-3\">
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/spider.pdf\" class=\"btn btn-primary btn-large\">Download PDF</a></p>
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/spider.png\" class=\"btn btn-primary btn-large\">Download PNG</a></p>
                                    </div>
                                    </div>
				    
				    <div class=\"col-md-11\">
                                    <div class=\"col-md-8\">
                                        <img class=\"img-thumbnail\" data-src=\"holder.js/500x500/auto\" alt=\"500x500\" src=\"circos.png\"><br>
                                        <p class=\"bg-info\">Genome-wide interaction circos map</p>
                                        <br><br>
                                    </div>
                                    <div class=\"col-md-3\">
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/circos.pdf\" class=\"btn btn-primary btn-large\">Download PDF</a></p>
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/circos.png\" class=\"btn btn-primary btn-large\">Download PNG</a></p>
                                    </div>
                                    </div>
				    
				    <div class=\"col-md-11\">
                                    <div class=\"col-md-8\">
                                        <img class=\"img-thumbnail\" data-src=\"holder.js/500x500/auto\" alt=\"500x500\" src=\"distance.png\"><br>
                                        <p class=\"bg-info\">Distance distribution to key features</p>
                                        <br><br>
                                    </div>
                                    <div class=\"col-md-3\">
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/distance.pdf\" class=\"btn btn-primary btn-large\">Download PDF</a></p>
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/distance.png\" class=\"btn btn-primary btn-large\">Download PNG</a></p>
                                    </div>
                                    </div>
                                    
				    <div class=\"col-md-11\">
                                    <div class=\"col-md-8\">
                                        <img class=\"img-thumbnail\" data-src=\"holder.js/500x500/auto\" alt=\"500x500\" src=\"DNA_replication.png\"><br>
                                        <p class=\"bg-info\">DNA replication timing of 4C regions</p>
                                        <br><br>
                                    </div>
                                    <div class=\"col-md-3\">
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/DNA_replication.pdf\" class=\"btn btn-primary btn-large\">Download PDF</a></p>
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/DNA_replication.png\" class=\"btn btn-primary btn-large\">Download PNG</a></p>
                                    </div>
				    </div>
                                    ";
        if (-e "/var/www/html/w4cseq/html/done/$id/$password/ChIP-Seq.pdf") {print HTML "
                                    <div class=\"col-md-11\">
				    <div class=\"col-md-8\">
                                        <img class=\"img-thumbnail\" data-src=\"holder.js/500x500/auto\" alt=\"500x500\" src=\"ChIP-Seq.png\"><br>
                                        <p class=\"bg-info\">Enrichment of user specified features</p>
                                        <br><br>
                                    </div>
                                    <div class=\"col-md-3\">
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/ChIP-Seq.pdf\" class=\"btn btn-primary btn-large\">Download PDF</a></p>
                                        <p><a href=\"http://w4cseq.usc.edu/done/$id/$password/ChIP-Seq.png\" class=\"btn btn-primary btn-large\">Download PNG</a></p>
                                    </div>
				    </div>
                                    ";}
        print HTML "
                                </div>
                                
                                <div class=\"col-md-9\" id = \"files\">
                                    <h2>Files</h2>
				    <p>Click to download</p>
                                    <div class=\"list-group\">
                                        <a href=\"http://w4cseq.usc.edu/done/$id/$password/FASTQ_FOR_MAPPING.fq\" class=\"list-group-item\"><strong>Fastq file trimmed for 4C analysis (For Illumina1.3+, base quality scores automatically converted to Sanger scores)</strong></a>
                                        <a href=\"http://w4cseq.usc.edu/done/$id/$password/FASTQ_FILTERED.fq\" class=\"list-group-item\"><strong>Fastq file filtered (mean base quality score for each read >=20) for alignment</strong></a>
                                        <a href=\"http://w4cseq.usc.edu/done/$id/$password/MAPPED_BAM.bam\" class=\"list-group-item\"><strong>Bam file of mapped reads</strong></a>
                                        <a href=\"http://w4cseq.usc.edu/done/$id/$password/MAPPED_BAM.bam.bai\" class=\"list-group-item\"><strong>Bam index file (you need this file to visualize bam result in IGV browser)</strong></a>
                                        <a href=\"http://w4cseq.usc.edu/done/$id/$password/DISTAL_INTERACTION_SITES.bed\" class=\"list-group-item\"><strong>Bed file of non-randomly mapped interacting sites with coverage</strong></a>
                                        <a href=\"http://w4cseq.usc.edu/done/$id/$password/positive_hits.bed\" class=\"list-group-item\"><strong>Bed file of significant interacting sites</strong></a>
					<a href=\"http://w4cseq.usc.edu/done/$id/$password/SIGNIFICANT_REGIONS.bed\" class=\"list-group-item\"><strong>Bed file of significant interacting regions</strong></a>
                                    </div>
                                </div>  
                            </article>
                                
                        </section>
                         
			<div class=\"col-md-12\" id = \"UCSC\">
			    <h2>View in UCSC genome browser</h2>
			    <span class=\"glyphicon glyphicon-log-in\"></span>
			    <a href=\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=$info{ref}&hgt.customText=http://w4cseq.usc.edu/done/$id/$password/UCSC_view.bed\">view in UCSC genome browser</a>
			</div>
			
			<br><br> 
                        <br><br>	
			<br><br>
			
                        <h4>
                           Please send questions or comments to <a href=\"mailto:caim\@usc.edu\">caim\@usc.edu</a>
                        </h4>
                        <hr>
                        <footer>
                           <p class=\"pull-right\"><a href=\"#\">Back to top</a></p>
                           <p>@ Wang genomics lab 2014</p>
                        </footer>
                        <br><br><br>
                    </div>
             
                    <script src=\"../html/W4CSEQ_files/js/jquery.min.js\" type=\"text/javascript\"></script>
                    <script src=\"../html/W4CSEQ_files/js/bootstrap.js\" type=\"text/javascript\"></script>
                </body>
	</html>
	";
	close (HTML);
	
	
	$email_header = "Dear W4CSEQ user, your submission (identifier: $id) was received at $info{submission_time} and processed at $process_time.\n";
	$email_header =~ s/(.{1,69})\s/$1\n/g;

	if ($failed_command) {
		$email_body = "We were unable to generate results for your submission due to an '$failed_command' error.\n";
	} else {
		$email_body = "Your submission is done: http://w4cseq.usc.edu/done/$id/$password/index.html\n\n";#### url
		$email_body .= "fastqfile=$info{query1}\nbuildver=$info{ref}\n\n";
		
		
		
		
	
		$email_tail .= "The citation for the above result is: http://w4cseq.usc.edu\n\n";
		$email_tail .= "Questions or comments may be directed to $CARETAKER.\n";
		$email_tail =~ s/(.{1,69})\s/$1\n/g;
	}
		
	open (EMAIL, ">$WORK_DIRECTORY/$id/email") or warn ">" . scalar (localtime) . " (id: $id)\ncannot create new mail $ARGV[1]\n" and return;
	flock (EMAIL, 2);
	print EMAIL "From: $CARETAKER\nReply-To: $CARETAKER\nSubject: 4CSEQ web server results for your query (identifier: $id)\n\n";
	print EMAIL $email_header, '-'x70, "\n\n", $email_body, '-'x70, "\n\n", $email_tail, "\n";
	flock (EMAIL, 2);
	close (EMAIL);

	if ($info{email}) {
		system ("/usr/sbin/sendmail $info{email} < $WORK_DIRECTORY/$id/email") and warn '>' . scalar(localtime) . " (id: $id)\ncannot send mail\n" and return 0;
	}
	
	#system("cp index.html $HTML_DIRECTORY/done/$id/$password") and warn "Error runnning <cp index.html $HTML_DIRECTORY/done/$id/$password>";
}


sub REAPER {
	my $pid;
	while (($pid = waitpid(-1, WNOHANG)) > 0) {
		$zombies--;
	}
}
