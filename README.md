clistats
================================================================================

clistats is a command line interface tool to compute statistics of a set
of delimited input numbers from a stream such as a Comma Separated Value (.csv)
or Tab Separated Value (.tsv) file. The default delimiter is a comma. Input data can
be a file, a redirected pipe, or the manually entered at the console. To stop
processing and display the statistics during manual input, enter the EOF signal
(CTRL-D on POSIX systems like Linux or Cygwin or CTRL-Z on Windows).

### I/O options

  * Input data can be from a file, standard input, or a pipe
  * Output can be written to a file, standard output, or a pipe
  * Output uses headers that start with "#" to enable piping to gnuplot
    
### Parsing options

  * Signal, end-of-file, or blank line based detection to stop processing
  * Comment and delimiter character can be set
  * Columns can be filtered out from processing
  * Rows can be filtered out from processing based on numeric constraint
  * Rows can be filtered out from processing based on string constraint
  * Rows can be sampled uniformly or randomly.
  * Initial header rows can be skipped
  * Fixed number of rows can be processed
  * Duplicate delimiters can be ignored
  * Rows can be reshaped into columns
  * Strictly enforce that only rows of the same size are processed
  * A row containing column titles can be used to title output statistics
    
### Statistics options

  * Summary statistics (Count, Minimum, Mean, Maximum, Standard deviation)
  * Covariance
  * Correlation
  * Least squares offset
  * Least squares slope
  * Histogram
  * Raw data after filtering
    
### Warnings
 - Delimiters are not preserved if in quotes. Any delimiter character will
   cause the row to be split.

 - All statistics are computed on a rolling basis and not by using the
   entire dataset, so all statistics (except minimum, maximum, and count) are
   approximations. When row sampling is used all statistics are approximations.

 - The histogram is also computed on a rolling basis with a
   dynamic histogram merging algorithm. New data points that do not fall
   within the current histogram's bounds are added to a cache. Once the cache
   is full it is merged with the current histogram. Bins sizes are scaled by
   the smallest integer needed to include the new cache data and maintain the
   same number of bins. The initial histogram is empty so all data is initially
   added to the cache. Therefore, the cache size should not be too small and
   preferably be set to a size that captures the statistics of the underlying
   sample. However, the larger the cache size the more memory required to
   store the cache values. The default value of the cache size is 1000.
   
### Alternatives
I trolled around online and searched for existing solutions to the same problem.
There are several very nice solutions which I've listed below. I still think
clistats is the most flexible, robust, and easiest to use out-of-the-box but
I've done no real testing so it's just coder's pride saying that. However,
there are some others that appear might be faster but are more limited in the
scope of what kind of input they can process and what output statistics they
generate.

 - [|STAT](http://hcibib.org/perlman/stat/)
 - [Average](http://sourceforge.net/projects/average/)
 - [datastat] (http://sourceforge.net/projects/datastat/)
 - [qstat](https://github.com/tonyfischetti/qstats)
 - [st] (https://github.com/nferraz/st)
 - [sta] (https://github.com/simonccarter/sta)
 - [stats](http://web.cs.wpi.edu/~claypool/misc/stats/stats.html)
 - [stats-tools] (https://github.com/jweslley/stats-tools)

EXAMPLES
================================================================================

Since there's no better explanation to running a tool like some simple examples,
below are some basic use cases to get you going on running clistats. The
following show how to provide input data, filter the data, and redirect the
computed statistics to gnuplot.

### Standard Input

Input data is taken from the standard input so a user can input
numbers at the console:

    ./clistats
    1,2,3,4
    5,6,7,8
    9,0,1,2
    3,4,5,6
    
    #============================================================================
    #                            Statistics
    #============================================================================
    #     Dimension     Count      Minimum         Mean      Maximum        Stdev
    #----------------------------------------------------------------------------
                  1         4     1.000000     4.500000     9.000000     2.958040
                  2         4     0.000000     3.000000     6.000000     2.236068
                  3         4     1.000000     4.000000     7.000000     2.236068
                  4         4     2.000000     5.000000     8.000000     2.236068

### File Input

An input file can be redirected to process delimited data from a file or by
specifying the input file (see -i option):

    ./clistats < file.csv

### Pipe Input

Input data can also be provided using a pipe:

    (echo "1,2,3,4"; echo "5,6,7,8"; echo "9,0,1,2"; echo "3,4,5,6") | ./clistats
    #============================================================================
    #                            Statistics
    #============================================================================
    #     Dimension     Count      Minimum         Mean      Maximum        Stdev
    #----------------------------------------------------------------------------
                  1         4     1.000000     4.500000     9.000000     2.958040
                  2         4     0.000000     3.000000     6.000000     2.236068
                  3         4     1.000000     4.000000     7.000000     2.236068
                  4         4     2.000000     5.000000     8.000000     2.236068

### Realistic Example

A slightly more realistic example would be to download some actual data. The
example below downloads comma delimited raw data from the [Lahman Baseball
Archive](http://www.seanlahman.com/baseball-archive/statistics/).
The results show various batting statistics over the years 1871 to 2013.

Columns that have entirely non-numeric data are displayed with a
Not-A-Number string "nan". The example makes use of displaying the data's
correlation table which can highlight trends between variables. For
instance, note the very low correlation of 0.253640 between Home Runs (HR)
and Stolen Bases (SB) as most power hitter don't tend to be faster runners.
Same goes for Triples (3B) as it's hard for those big guys to make it all
the way to 3rd base so their correlation is 0.338364.

Note, you may need to install wget and unzip to get this example to work.
Alternatively, you can download and unzip the files manually. Also, I don't
own any of this archive data so use it within the site's legal provisions.

    wget http://seanlahman.com/files/database/lahman-csv_2014-02-14.zip
    unzip lahman-csv_2014-02-14.zip
    ./clistats --titles 1 --filterColumn "1,8:10,12:13,15" --correlation < Batting.csv
    #=============================================================================
    #                           Correlation
    #=============================================================================
    #   playerID         AB          R          H         3B         HR         SB
    #-----------------------------------------------------------------------------
             nan        nan        nan        nan        nan        nan        nan
             nan   1.000000   0.950196   0.987135   0.712319   0.684625   0.603282
             nan   0.950196   1.000000   0.965945   0.742781   0.719900   0.657723
             nan   0.987135   0.965945   1.000000   0.736148   0.693786   0.611282
             nan   0.712319   0.742781   0.736148   1.000000   0.338364   0.609333
             nan   0.684625   0.719900   0.693786   0.338364   1.000000   0.253640
             nan   0.603282   0.657723   0.611282   0.609333   0.253640   1.000000

### Plotting with gnuplot

You can pipe this output to gnuplot:

    ./clistats --titles 1 --filterColumn "1,8:10,12:13,15" --correlation < Batting.csv | gnuplot -p -e 'plot "-" matrix with image title "Correlation"'

which will display an image representation of the correlation matrix. You'll
see any "nan" strings are interpreted as undefined values by gnuplot and
rendered as black.

![Correlation](/examples/correlation.png?raw=true "Correlation")

### Column filtering

Column filters can be used to removed unwanted columns from the computed
statistics using the "--filterColumn" option. The example above can be
modified to additionally filter out the first column and remove the
unwanted "nan" strings.

    ./clistats --titles 1 --filterColumn "8:10,12:13,15" --correlation < Batting.csv | gnuplot -p -e 'plot "-" matrix with image title "Correlation"'

![Column Filtering](/examples/column-filtering.png?raw=true "Column Filtering")

### Row filtering

Row filters can also be used to remove unwanted rows from the computed statistics
using either numeric or string criteria. Multiple row filters can be used and will
be processed in the order provided. However, all string filters will be processed
first then all numeric filters are processed. The following example keeps only
those batting statistics after the year 2000 by matching the entries in the 2nd
column to the interval [2000,infinity]:

    ./clistats --titles 1 --filterColumn "2,8:10,12:13,15" --filterNumeric "2,2000,inf" < Batting.csv
    #===========================================================================
    #                           Statistics
    #===========================================================================
    #   Dimension   Count       Minimum          Mean       Maximum        Stdev
    #---------------------------------------------------------------------------
           yearID   18641   2000.000000   2006.440051   2013.000000     4.028313
               AB   17377      0.000000    134.062611    716.000000   186.083370
                R   17377      0.000000     18.102664    152.000000    28.139381
                H   17377      0.000000     35.189561    262.000000    52.407670
               3B   17377      0.000000      0.731369     23.000000     1.657640
               HR   17377      0.000000      4.080566     73.000000     7.865398
               SB   17377      0.000000      2.308741     78.000000     6.119866

The following example adds an additional filter to keep only those batting
statistics for players from Boston by matching the string "BOS" to the 4th column:

    ./clistats --titles 1 --filterColumn "2,8:10,12:13,15" --filterNumeric "2,2000,inf" --filterString "4,BOS" < Batting.csv
    #===========================================================================
    #                           Statistics
    #===========================================================================
    #   Dimension   Count       Minimum          Mean       Maximum        Stdev
    #---------------------------------------------------------------------------
           yearID     644   2000.000000   2006.335404   2013.000000     4.014638
               AB     543      0.000000    145.392265    660.000000   197.596551
                R     543      0.000000     21.965009    123.000000    32.242187
                H     543      0.000000     39.930018    213.000000    57.410212
               3B     543      0.000000      0.720074     13.000000     1.586924
               HR     543      0.000000      4.974217     54.000000     8.855879
               SB     543      0.000000      2.123389     70.000000     6.290502

USAGE
================================================================================

To display a full listing of the application options use

    $ ./clistats --help

INSTALL
================================================================================

###Build with Make

A simple GNU make file is provided to build the code. Just run "make".

###Build with CMake

A CMake file is also provided to build out of source and provided for
future development.

  * Definitions
    1. \<source>    Directory where source code was installed
    2. \<build>     Directory where code will be built
    3. \<install>   Directory where executable will be installed
  * Prerequisites:
    1. CMake 2.8 (or higher)            To build from the Cmake files.
    2. Visual Studio 2008 (or higher)   To build on Windows
  * Instructions:
    1. Create build directory: mkdir \<build>
    2. Change to build directory: cd \<build>
    3. Build
      - Linux
        * cmake \<source> -DCMAKE_INSTALL_PREFIX=\<install> -DCMAKE_BUILD_TYPE=Release
        * make && make install
      - Windows
        * cmake \<source> -DCMAKE_INSTALL_PREFIX=\<install>
        * Open Visual Studio solution "clistats.sln"
        * Run project "ALL_BUILD"
        * Run project "INSTALL"

LICENSE
================================================================================

Copyright (c) 2014 Daniel Pulido <dpmcmlxxvi@gmail.com>

clistats is released under the [MIT License](http://opensource.org/licenses/MIT)

CHANGELOG
================================================================================

- Version 0.1
    
  * Initial release

AUTHOR
================================================================================

Copyright 2014 by Daniel Pulido <dpmcmlxxvi@gmail.com>
