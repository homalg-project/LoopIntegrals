# SPDX-License-Identifier: GPL-2.0-or-later
# LoopIntegrals: Compute master integrals using commutative and noncommutative methods from computational algebraic geometry
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "LoopIntegrals",
Subtitle := "Compute master integrals using commutative and noncommutative methods from computational algebraic geometry",
Version := Maximum( [
                   "2022.04-01", ## Mohamed's version
                   ## this line prevents merge conflicts
                   ] ),

Date := "26/04/2022",

License := "GPL-2.0-or-later",


Persons := [
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Mohamed",
    LastName := "Barakat",
    WWWHome := "https://mohamed-barakat.github.io/",
    Email := "mohamed.barakat@uni-siegen.de",
    PostalAddress := Concatenation(
               "Walter-Flex-Str. 3\n",
               "57068 Siegen\n",
               "Germany" ),
    Place := "Siegen",
    Institution := "University of Siegen",
  ),
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Robin",
    LastName := "Brüser",
    WWWHome := "https://www.researchgate.net/scientific-contributions/2138634848_Robin_Brueser",
    Email := "brueser@uni-mainz.de",
    PostalAddress := Concatenation(
               "Walter-Flex-Str. 3\n",
               "57068 Siegen\n",
               "Germany" ),
    Place := "Siegen",
    Institution := "University of Siegen",
  ),
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Tobias",
    LastName := "Huber",
    WWWHome := "https://www.physik.uni-siegen.de/tp1/research/researchgroups/huber.html",
    Email := "huber@tp1.physik.uni-siegen.de",
    PostalAddress := Concatenation(
               "Walter-Flex-Str. 3\n",
               "57068 Siegen\n",
               "Germany" ),
    Place := "Siegen",
    Institution := "University of Siegen",
  ),
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Jan",
    LastName := "Piclum",
    WWWHome := "https://www.researchgate.net/scientific-contributions/34408173_Jan_Piclum",
    Email := "piclum@physik.uni-siegen.de",
    PostalAddress := Concatenation(
               "Walter-Flex-Str. 3\n",
               "57068 Siegen\n",
               "Germany" ),
    Place := "Siegen",
    Institution := "University of Siegen",
  ),
],

# BEGIN URLS
SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/homalg-project/LoopIntegrals",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome  := "https://homalg-project.github.io/pkg/LoopIntegrals",
PackageInfoURL  := "https://homalg-project.github.io/LoopIntegrals/PackageInfo.g",
README_URL      := "https://homalg-project.github.io/LoopIntegrals/README.md",
ArchiveURL      := Concatenation( "https://github.com/homalg-project/LoopIntegrals/releases/download/v", ~.Version, "/LoopIntegrals-", ~.Version ),
# END URLS

ArchiveFormats := ".tar.gz .zip",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=  "",

PackageDoc := rec(
  BookName  := "LoopIntegrals",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Compute master integrals using commutative and noncommutative methods from computational algebraic geometry",
),

Dependencies := rec(
  GAP := ">= 4.9.1",
  NeededOtherPackages := [
                   [ "GAPDoc", ">= 1.5" ],
                   [ "MatricesForHomalg", ">= 2020.05.20" ],
                   [ "RingsForHomalg", ">= 2020.10-05" ],
                   [ "GradedRingForHomalg", ">= 2020.05.01" ],
                   ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := function()
        return true;
    end,

TestFile := "tst/testall.g",

Keywords := [ "loop diagrams", "generating vectors", "syzygies", "master integrals" ],

));
