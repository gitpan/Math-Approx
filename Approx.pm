#!/usr/local/ls6/perl/bin/perl
#                              -*- Mode: Perl -*- 
# Approx.pm -- 
# ITIID           : $ITI$ $Header $__Header$
# Author          : Ulrich Pfeifer
# Created On      : Wed Oct 25 10:50:38 1995
# Last Modified By: Ulrich Pfeifer
# Last Modified On: Wed Oct 25 12:47:41 1995
# Language        : Perl
# Update Count    : 57
# Status          : Unknown, Use with caution!
# 
# (C) Copyright 1995, Universität Dortmund, all rights reserved.
# 
# $Locker:  $
# $Log: Approx.pm,v $
# Revision 1.2  1995/10/25  12:38:59  pfeifer
# Added documentation.
#
# Revision 1.1  1995/10/25  10:29:26  pfeifer
# Initial revision
#
# 

=head1 NAME

Math::Approx

=head1 METHODS

=head2 new

    new Math::Approx (\&poly, 5, %x);

The first argument after the class name must be a reference to
function which takes two arguments: The I<degree> and the I<x> value.

For interpolation with plain polynomials I<poly> can be defined as:

        sub poly {
            my($n,$x) = @_;
        
            return $x ** $n;
        }

The second argument is the maximum degree which should be used for
interpolation. Degrees start with B<0>. 

The rest of the arguments are treated as pairs of B<x> and B<y>
samples which should be approximated.

The method returns a Math::Approx reference.

=head2 approx

	$approximation->approx(17);

The method returns the approximated  B<y> value for the B<x> value
given as argument.

=head2 fit

 	$approximation->fit;

Returns the medim square error for the data points.

=head2 plot

 	$approximation->plot("tmp/app");

Prints all data pairs and the corresponding approximation pairs in a
file whichs filename is given as argument. The file should be suitable
for usage with gnuplot(1).


=head2 print

 	$approximation->print;

Prints information about the approximation on I<STDOUT>

=head1 EXAMPLE

        use Math::Approx;
        
        sub poly {
            my($n,$x) = @_;
        
            return $x ** $n;
        }
        
        for (1..20) {
            $x{$_} = sin($_/10)*cos($_/30)+0.3*rand;
        }
        
        $a = new Math::Approx (\&poly, 5, %x);
        $a->print;
        $a->plot("mist");
        print "Fit: ", $a->fit, "\n";

=head1 SEE ALSO

gnuplot(1).

=head1 AUTHOR

Ulrich Pfeifer <pfeifer@ls6.informatik.uni-dortmund.de>

=cut


package Math::Approx;
use Math::Matrix;

$RCS_Id = '$Id: Approx.pm,v 1.2 1995/10/25 12:38:59 pfeifer Exp $ ';
($my_name, $my_version) = $RCS_Id =~ /: (.+).pm,v ([\d.]+)/;

sub new {
    my $type = shift;
    my $func = shift;
    my $degr = shift;
    my %data = @_;
    my $self = {};
    my @m;
    my @x;

    $self->{'F'} = $func;       # function of two arguments
    $self->{'N'} = $degr;
    $self->{'D'} = \%data;

    for $x (keys %data) {
        my $row = [];
        for $n (0 .. $degr) {
            push @{$row}, &{$func}($n, $x);
        }
        push @x, $data{$x};
        push(@m, $row);
    }
    my $I = new Math::Matrix(@m);            # $I->print("Initial\n");
    my $T = $I->transpose->multiply($I);     # $T->print("Quadratic\n");
    my $x = new Math::Matrix(\@x);           # $x->print("X\n");
    my $tx = $I->transpose->
        multiply($x->transpose);             # $tx->print("TX\n");
    my $eq = $T->concat($tx);                # $eq->print("EQN\n");
    my $s = $eq->solve;                      # $s->print("SOLUTION\n");
                                             # $T->multiply($s)->print("TEST\n");
    $self->{'A'} = $s->transpose->[0];
    bless $self, $type;
}

sub approx {
    my $self = shift;
    my $x    = shift;
    my $func = $self->{'F'};
    my $degr = $self->{'N'};
    my $result;

    for $n (0 .. $degr) {
        $result += &{$func}($n, $x) * $self->{'A'}->[$n];
    }

    $result;
}

sub fit {
    my $self = shift;
    my $result;
    my $n;

    for $key (keys %{$self->{'D'}}) {
        $result += ($self->{'D'}->{$key}-$self->approx($key))**2;
        #print STDERR "## $result\n";
        $n++;
    }

    $result/$n;
}

sub print {
    my $self = shift;
    
    printf "Function: %s\n", $self->{'F'};
    printf "Degree:   %d\n", $self->{'N'};
    print  "Koeff: ", join(' ', @{$self->{'A'}}), "\n";
    print  "Fit: ",    $self->fit, "\n";
    print  "Data:\n";
    print  "     X          Y     Approximation\n";
    for $key (sort {$a <=> $b} keys %{$self->{'D'}}) {
        printf("%10.5f %10.5f %10.5f\n", $key, $self->{'D'}->{$key}, 
               $self->approx($key));
    }
}

sub plot {
    my $self = shift;
    my $file = shift;
    open(OUT, ">$file") || die "Could not open $file: $!\n";
    
    print OUT "\n#data\n";
    for $key (sort {$a <=> $b} keys %{$self->{'D'}}) {
        print OUT $key, ' ', $self->{'D'}->{$key}, "\n";
    }

    print OUT "\n#Approximation\n";
    for $key (sort {$a <=> $b} keys %{$self->{'D'}}) {
        print OUT $key, ' ', $self->approx($key), "\n";
    }
    close(OUT);
    1;
}

1;
