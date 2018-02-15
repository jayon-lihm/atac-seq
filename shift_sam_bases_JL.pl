#!/usr/bin/perl

$plus_shift = 4;
$minus_shift = -5;


$chrfile = $chrfile = "/seq/jlihm/resources/mm9/mm9.chrom.sizes";
#$infile = $ARGV[0];
#$outfile = $ARGV[1];

#if ($#ARGV ne 1) {
#  print "commond line: perl shift_sam_bases.pl input_file output_file\n";
#  exit;
#}

open (in, "<$chrfile");
while ($line=<in>) {
  chomp $line;
  @data = split /\t/, $line;
  $chr = $data[0];
  $chrsize{$chr} = $data[1];
}
close in;

#open (in, "<$infile");
#open (out, ">$outfile");

while ($line=<>) {
 
  $row_count ++;
  #if ($row_count%1000000 eq 0) {
  #  print "$row_count\n";
  #}

  chomp $line;
  @data = split /\t/, $line;

  if ( substr($line, 0, 1) ne '@' ){
      $chr = $data[2];
      $plus_posi = $data[3];
      $minus_posi = $data[7];
      $length = $data[8];
      
      $flag = $data[1];
      $flagbin = dec2bin($flag);
      @flagbin = split "", $flagbin;
      
      ## Set $flagbin to be at least 5 digits in binary
      $flagbin_len = length $flagbin;
      if ($flagbin_len < 5){
	  $fillin = 5 - $flagbin_len;
	  $fillin_char ='0' x $fillin;
	  $flagbin = "$fillin_char$flagbin";
      }
      
      ##FLAG: 4 0x4 segment unmapped
      ##binary: 0000100
      ## "$flagbin[-3]==0" ==> replaced with "substr"
      ##FLAG: 16 0x10 SEQ being reverse complemented
      ##binary: 00010000
      ## "$flagbin[-5]==0" ==> replaced with "substr"
      $flagbin_l3 = substr($flagbin, -3, -2);
      $flagbin_l5 = substr($flagbin, -5, -4);
      
      if (($flagbin_l3 eq 0) and ($flagbin_l5 eq 0)) { 
	  $data[3] = $plus_posi + $plus_shift;
      }
      elsif (($flagbin_l3 eq 0) and ($flagbin_l5 eq 1)) {
	  $data[3] = $plus_posi + $minus_shift;
      }
  
      if (($data[3] > 0) and ($data[3] < $chrsize{$chr})) {
	  $pline = join "\t", @data;
	  print "$pline\n";
      }
  } else {
      print "$line\n";
  }
}

close in;
close out;


$decimalnum = bin2dec($binarynum);
# Convert a binary number to a decimal number
sub bin2dec {
  unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}


$binarynum = dec2bin($decimalnum);
# Convert a decimal to a binary
sub dec2bin {
  my $str = unpack("B32", pack("N", shift));
  $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
  return $str;
}

