#                              -*- Mode: Perl -*- 
# $Basename: Makefile.PL $
# $Revision: 1.2 $
# Author          : Ulrich Pfeifer
# Created On      : Fri Oct 27 11:19:11 1995
# Last Modified By: Ulrich Pfeifer
# Last Modified On: Tue Mar  3 15:15:40 1998
# Language        : Perl
# 
# (C) Copyright 1995,1998, Ulrich Pfeifer
# 


use ExtUtils::MakeMaker;

&WriteMakefile
  (
   'NAME'          => 'Math::Approx',
# $Format: "    'VERSION' => sprintf('%5.3f', $ProjectVersion$),"$
    'VERSION' => sprintf('%5.3f', 0.2),
   'dist'          => {
                       COMPRESS => 'gzip -9f',
                       PREOP    => 'pod2text Approx.pm > README',
                      },
   'clean' 	   => {'FILES' => '*~'},
   'PREREQ_PM'     => { 'Math::Matrix' => 0},
  );
