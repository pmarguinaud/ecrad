#!/usr/bin/perl -w


use strict;

while (<>)
  {
    next if (m/^(?:\s*\w+\s*=)?\s*(?: 0,)*\s*(?: 0 ;)?\s*$/o);
    next if (m/^\s*float\s+.*;$/o);
    print;
  }
