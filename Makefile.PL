#!/usr/local/ls6/perl/bin/perl
#                              -*- Mode: Perl -*- 
# Makefile.PL -- 
# ITIID           : $ITI$ $Header $__Header$
# Author          : Ulrich Pfeifer
# Created On      : Fri Oct 27 11:19:11 1995
# Last Modified By: Ulrich Pfeifer
# Last Modified On: Fri Oct 27 11:39:49 1995
# Language        : Perl
# Update Count    : 8
# Status          : Unknown, Use with caution!
# 
# (C) Copyright 1995, Universität Dortmund, all rights reserved.
# 
# $Locker$
# $Log$
# 


use ExtUtils::MakeMaker;
use Config;

$Version = '0.1';

&WriteMakefile(
               'NAME'     => 'Math::Approx',
               'LINKTYPE' => '$(INST_PM)',,
               'VERSION'  => $Version,
               'dist'     => { SUFFIX => 'gz', COMPRESS => 'gzip -f' },
               'clean' 	  => {'FILES' => '*~'},
);
