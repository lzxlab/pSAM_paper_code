#!/usr/bin/perl
### Author:    Lin Sheng, Ph.D.
###  Group:    SZIPP, Shenzhen, China, 518040
###  Email:    ls03ls@mail2.sysu.edu.cn

#-------------------------------------module--------------------------------------------------------

use warnings;
use strict;

#-------------------------------------global--------------------------------------------------------

my $func_name="sgRNA Extracter Forward Strands";
my $create_date="2017-05-04";
my $revise_date="2017-08-15";
my $program_ver="2.0";

my $str1="AAAGGACGAAACACCG";
my $str2="GTTTTAGAGCTAGAAATAGC";
my $str3="";
my $str4="";
my $infile="test_in.fq";
my $outfile="test_out.fq";

#-------------------------------------mainbody------------------------------------------------------

my ($l,@s,%h);
&show_welcome;
&get_io;
&run;
&quit;

#--------------------------------------subrountines-------------------------------------------------

sub show_welcome
{
    system "cls";
    my $sep="="x72;
    print "\n$sep\n\n";
    print "\t\t< $func_name >\n\n" if ($func_name ne "");
    print "\tAuthor_Name:       Lin Sheng Ph.D.\n";
    print "\tAffiliation:       SZIPP, Shenzhen, China\n";
    print "\tCreate_Date:       $create_date\n" if ($create_date ne "");
    print "\tRevise_Date:       $revise_date\n" if ($revise_date ne "");
    print "\tProgram_Ver:       $program_ver\n" if ($program_ver ne "");
    print "\n$sep\n\n";
}

sub status
{
    my $string=$_[0];
    my $status="[$_[1]]";
    printf "Status: %-50s%-10s\n",$string,$status;
}

sub error
{
    my $string=$_[0];
    printf "Errors: %s\n>>> Press enter key to exit..",$string;
    exit;
}

sub quit
{
    my $sec=2;
    print "Status: mission success! Waiting $sec seconds to Exit..\n";
    sleep($sec);
    exit;
}

#--------------------------------------subrountines-------------------------------------------------

sub get_io
{
    print ">>> Pls enter infile name (eg. test_in.fq):\n";
    $infile=$ARGV[0];
    chomp $infile;
    print ">>> Pls enter outfile name (eg. test_out.fq):\n";
    $outfile=$ARGV[1];
    chomp $outfile;
}

sub get_rev
{
    my $str=shift;
    my @str=split(//,$str);
    my @rev=reverse(@str);
    my $newstr="";
    foreach my $item (@rev)
    {
        my $com;
        if($item eq "A")
        {
            $com="T";
        }
        elsif($item eq "T")
        {
            $com="A";
        }
        elsif($item eq "G")
        {
            $com="C";
        }
        elsif($item eq "C")
        {
            $com="G"
        }
        else
        {
            $com="N";
        }
        $newstr.=$com;
    }
    return $newstr;
}

sub just_rev
{
    my $str=shift;
    my @str=split(//,$str);
    my @rev=reverse(@str);
    my $newstr=join("",@rev);
    return $newstr;
}

sub run
{
    open (IN,"$infile") or die "cannot open infile!";
    open (OUT,">$outfile") or die "cannot open outfile!";
    
    $str3=&get_rev($str2);
    $str4=&get_rev($str1);
    
    # dealing with each record
    my $n1=0;
    my $n2=0;
    while($l=<IN>)
    {
        my $symbol=$l;
        my $seq=<IN>;
        my $strand=<IN>;
        my $qual=<IN>;
        if($seq=~m/$str1(\w+)$str2/)
        {
            $n1++;
            my $tarseq=$1;
            my $pos=index($seq,$tarseq);
            my $len=length($tarseq);
            my $tarqual=substr($qual,$pos,$len);
            print OUT "$symbol";
            print OUT "$tarseq\n";
            print OUT "$strand";
            print OUT "$tarqual\n";
        }
        elsif($seq=~m/$str3(\w+)$str4/)
        {
            $n2++;
            my $tarseq=$1;
            my $pos=index($seq,$tarseq);
            my $len=length($tarseq);
            my $tarqual=substr($qual,$pos,$len);
            
            # trans to fwd strand
            $tarseq=&get_rev($tarseq);
            $tarqual=&just_rev($tarqual);
            $strand="-\n";
            print OUT "$symbol";
            print OUT "$tarseq\n";
            print OUT "$strand";
            print OUT "$tarqual\n";
        }
    }
    close (IN);
    close (OUT);
    &status("Found fwd records:",$n1);
    &status("Found rev records:",$n2);
}

































