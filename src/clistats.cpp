
/**
 * Computes command line interface statistics for a stream of delimited input numbers.
 * @author Daniel Pulido <dpmcmlxxvi@gmail.com>
 * @copyright Copyright (c) 2014 Daniel Pulido <dpmcmlxxvi@gmail.com>
 * @file clistats.cpp
 * @license MIT License (http://opensource.org/licenses/MIT)
 */

/* Standard headers */
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

// macros for logging message
#define LOG_MESSAGE(level, message)             \
    {                                           \
        if (level <= Logger::logLevel)          \
        {                                       \
            std::stringstream msg;              \
            msg << message << std::endl;        \
            Logger::log(level) << msg.str();    \
        }                                       \
    }

/**
 * @struct ApplicationProperties
 * @brief Application wide properties
 */
struct ApplicationProperties
{
    /**
     * Application version
     * @returns Version string
     */
    static
    std::string
    version()
    {
        return ApplicationProperties::VERSION_MAJOR + "." +
            ApplicationProperties::VERSION_MINOR + "." +
            ApplicationProperties::VERSION_PATCH;
    }
    static std::string NAME;
    static std::string AUTHOR;
    static std::string VERSION_MAJOR;
    static std::string VERSION_MINOR;
    static std::string VERSION_PATCH;
};

std::string ApplicationProperties::NAME = "clistats";
std::string ApplicationProperties::AUTHOR = "dpmcmlxxvi@gmail.com";
std::string ApplicationProperties::VERSION_MAJOR = "1";
std::string ApplicationProperties::VERSION_MINOR = "0";
std::string ApplicationProperties::VERSION_PATCH = "0";

/**
 * @struct Logger
 * @brief Simple logger
 * @details Should be invoked via macro "LOG_MESSAGE(Logger::Level, "message");"
 */
struct Logger
{

    /**
     * @enum Logger level
     */
    struct Level
    {
        /**
         * Level type
         */
        enum Type
        {
            FATAL = 0,
            ERROR,
            WARNING,
            INFO,
            DEBUG,
            DETAIL
        };
        
    };
    
    /**
     * Log message. Logger level is written to current stream line then returns stream. User is expected to add EOL.
     * @return Logging stream corresponding to log level
     */
    static std::ostream & log(Logger::Level::Type level)
    {
        std::string label = (level==Logger::Level::FATAL  ? "FATAL":
                            (level==Logger::Level::ERROR  ? "ERROR":
                            (level==Logger::Level::WARNING? "WARNING":
                            (level==Logger::Level::INFO   ? "INFO":
                            (level==Logger::Level::DEBUG  ? "DEBUG":"DETAIL")))));
        std::ostream * os = (level <= Logger::Level::WARNING ? &std::cerr : &std::cout);
        *os << std::left << std::setw(7) << label << ": ";
        return *os;
    }

    /**
     * Global application log level
     */
    static Level::Type logLevel;
    
};

// Initialize static variable
Logger::Level::Type Logger::logLevel = Logger::Level::FATAL;

/**
 * @class StringParser
 * @brief String parsing methods
 */
class StringParser
{

public:

    /**
     * Parse number to string
     * @param[in] value Value to parse
     * @param[in] isScientific True of string format should be scientific otherwise format is fixed
     * @return String representation of input value
     */
    template <class T>
    static
    std::string
    parseNumber(T value,
                bool isScientific = false)
    {
        // Check for NaN to consistently return "nan"
        volatile double d = (double)value;
        if (d != d) return "nan";
        std::string entry;
        std::stringstream parser;
        parser << (isScientific ? std::scientific : std::fixed) << value;
        parser >> entry;
        return entry;
    }

    /**
     * Parse statistic to string. If statistic was based on a count of zero then returns nan string.
     * @param[in] count Number of data points on which statistic was based on.
     * @param[in] value Value to parse
     * @param[in] nan String to return is statistic was not based on non-zero counts
     * @param[in] isScientific True of string format should be scientific otherwise format is fixed
     * @return String representation of input statistic
     */
    template <class T>
    static
    std::string
    parseStatistic(int count,
                   T value,
                   std::string nan = "nan",
                   bool isScientific = false)
    {
        if (count <= 0) return nan;
        return StringParser::parseNumber<T>(value, isScientific);
    }

    /**
     * Replace first n occurences of src with dst in content
     * @param[out] content String to replace content of
     * @param[in] src Source string to search for and replace with dst
     * @param[in] dst Destination string to replace src with
     * @param[in] n Number of times with which to replace src
     */
    static
    void
    replacen(std::string & content,
             std::string const & src,
             std::string const & dst,
             unsigned int const n = 1)
    {
        if(src.empty()) return;

        unsigned int count = 0;
        unsigned int pos = 0;
        while((pos = content.find(src, pos)) != std::string::npos)
        {
            content.replace(pos, src.length(), dst);
            pos += dst.length();
            if (++count == n) return;
        }
    }

    /**
     * Convert string to template type value
     * @param[in] src Source string
     * @param[out] dst Destination value
     * @return True if parsing successful otherwise false
     */
    template <class T>
    static
    bool
    toValue(std::string const & src,
            T & dst)
    {
        std::istringstream parser(src.c_str());
        parser >> dst;
        if (!parser) return false;
        return (parser.rdbuf()->in_avail() == 0);
    }

    /**
     * Trim leading and trailing white spaces
     * @param[in] src Source string
     * @return Trimmed string
     */
    static
    const std::string
    trim(const std::string & src)
    {
        std::string::size_type startpos = src.find_first_not_of(" \t");
        std::string::size_type endpos = src.find_last_not_of(" \t");
        return (std::string::npos==startpos ? "" : src.substr(startpos,endpos-startpos+1));
    }

    /**
     * Update width with current value's string size
     * @param[in] value Value to parse
     * @param[inout] width Width to update
     */
    template <class T>
    static
    void
    updateWidthNumber(T value,
                      long & width)
    {
        std::string entry = StringParser::parseNumber<T>(value);
        width = (width < (long) entry.size() ? entry.size() : width);
    }

    /**
     * Update width with current value's string size
     * @param[in] value Value to parse
     * @param[inout] width Width to update
     */
    static
    void
    updateWidthString(std::string value,
                      long & width)
    {
        std::string entry = value;
        width = (width < (long) entry.size() ? entry.size() : width);
    }

};

/**
 * @class StringSplitter
 * @brief Splits a source string into tokens given a delimiter character
 */
class StringSplitter
{

public:

    /**
     * Subscript type for use with operator[]
     */
    typedef std::vector<std::string>::size_type size_type;

public:

    /**
     * Create StringSplitter
     */
    StringSplitter()
    {
    }

    /**
     * Create and initialize a new StringSplitter
     * @param[in] src Source string to split
     * @param[in] delim The delimiter to split the string around
     * @param[in] duplicate If true then duplicate sequential delimiters are treated as one delimiter
     * @param[in] mask Boolean mask of tokens to keep. If token index is false in mask or is out of range then that token not discarded.
     */
    StringSplitter(std::string const & src,
                   std::string const & delim,
                   bool duplicate = false,
                   std::vector<bool> const & mask = std::vector<bool>())
    {
        this->initialize(src, delim, duplicate, mask);
    }

    /**
     * Retrieve the token at the specified index
     * @param[in] index Token index
     * @return The token at the specified index
     */
    std::string const &
    at(size_type index) const
    {
        return this->_tokens.at(index);
    }

    /**
     * Get reference to underlying split string vector
     * @return Vector of split string
     */
    const std::vector<std::string> &
    get() const
    {
        return this->_tokens;
    }
    
    /**
     * Initialize splitter with a new source and delimiter
     * @param[in] src The string to split
     * @param[in] delim The delimiter to split the string around
     * @param[in] duplicate If true then duplicate sequential delimiters are treated as one delimiter
     * @param[in] mask Boolean mask of tokens to keep. If token index is false in mask or is out of range then that token not discarded.
     */
    void
    initialize(std::string const & src,
               std::string const & delim,
               bool duplicate = false,
               std::vector<bool> const & mask = std::vector<bool>())
    {
        std::vector<std::string> tokens;
        std::string::size_type startpos = 0;
        std::string::size_type endpos = 0;
        unsigned int index = 0;
        bool doParse = true;
        while (doParse)
        {

            endpos = src.find(delim, startpos);

            // Check delimiter is not duplicate
            bool isDuplicate = (duplicate && (endpos == startpos));
            if (!isDuplicate)
            {
                // Check entry is not masked out
                bool isMasked = (mask.size()!=0) && (( (index >= mask.size()) || (!mask.at(index)) ));
                if (!isMasked)
                {
                    tokens.push_back(src.substr(startpos, endpos - startpos));
                }
                index++;
            }

            // We just copied the last token
            if (endpos == std::string::npos) break;

           // Exclude the delimiter in the next search
            startpos = endpos + delim.size();

        }
        this->_tokens.swap(tokens);
    }

    /**
     * Retrieve the number of split tokens
     * @return The number of split tokens
     */
    size_type
    size() const
    {
        return this->_tokens.size();
    }

    /**
     * Retrieve all the tokens
     * @return Vector of split tokens
     */
    std::vector<std::string>
    tokens() const
    {
        return this->_tokens;
    }

    /**
     * Convert string to template type value
     * @param[in] index Token index
     * @param[out] dst Destination value
     * @return True if parsing successful otherwise false
     */
    template <class T>
    bool
    toValue(size_type index,
            T & dst) const
    {
        return StringParser::toValue<T>(this->at(index), dst);
    }

    /**
     * Convert string to template type value
     * @param[in] index Token index
     * @param[out] dst Destination value
     * @return True if parsing successful otherwise false
     * @throws std::runtime_error if unable to parse string
     */
    template <class T>
    T
    toValue(size_type index) const
    {
        T & dst;
        bool isGood = this->toValue<T>(this->at(index), dst);
        if (!isGood) throw std::runtime_error("Unable to parse string = " + this->at(index));
    }
    
    /**
     * Convert string to range of integers. Expected format is
     * "A:B". If a range is found (e.g., 3:5) then the array
     * of numbers from the first to the second are returned
     * (e.g., [3,4,5]). If no range delimiter is found then the string
     * is assumed to be a single integer and parsing is attempted.
     * @param[in] src Source string
     * @param[out] dst Destination vector to which integers are appended
     * @param[in] delim Optional range delimiter string. Default = ":"
     * @return True if parsing succeeded otherwise false
     */
    static
    bool
    toIntegers(std::string const & src,
               std::vector<int> & dst,
               std::string delim = ":")
    {
    
        int start = 0;
        int stop = 0;

        StringSplitter splitter(src, delim);
            
        if (splitter.size() == 1)
        {
            // If one token then try to parse and set start = stop
            bool good = StringParser::toValue<int>(splitter.at(0), start);
            if (!good) return false;
            stop = start;
        }
        else if (splitter.size() == 2)
        {
            // If two tokens then try to parse start and stop
            bool goodStart = StringParser::toValue<int>(splitter.at(0), start);
            if (!goodStart) return false;
            
            bool goodStop = StringParser::toValue<int>(splitter.at(1), stop);
            if (!goodStop) return false;
            
            // Swap so start is before stop
            if (stop < start)
            {
                std::swap(start, stop);
            }
        }
        else            
        {
            // If not 1 or 2 tokens then fail
            return false;
        }

        // Create array of integers
        for (int i = start; i < stop+1; i++)
        {
            dst.push_back(i);
        }
        
        return true;
        
    }
    
private:

    /**
     * Contains the split tokens
     */
    std::vector<std::string> _tokens;

};

/**
 * @class FixedSizeCache
 * @brief Implements a fixed size cache of sortable items. Once full it is up to the user
 *        to empty the cache with #reset before adding more items.
 * @details As items are added the cache count increases. However, the cache count
 *          is different from its size. The size is fixed while the count will satisfy
 *          0 <= count <= size. In addition, the smallest and largest items are computed.
 *          Therefore, the comparison operator "<" must be defined for the item type stored.
 */
template <class T>
class FixedSizeCache
{

public:

    /**
     * @param[in] size Number of cache entries to store before merging with histogram. Must be larger than 0.
     */
    FixedSizeCache(int const size = 1000) :
        _count(0)
    {
        this->initialize(size);
    }

    /**
     * Add value to cache
     * @throws std::exception if cache is already full
     */
    void
    add(T value)
    {
        if (this->full()) throw std::runtime_error("FixedSizeCache::add() called on full cache");
        this->_data.at(this->_count) = value;
        this->_count++;
    }

    /**
     * Get cache entry
     * @throws std::exception if index if out-of-range of cache count
     * @return Cache entry at index
     */
    T
    at(int const index) const
    {
        if (index >= this->_count) throw std::runtime_error("FixedSizeCache index out-of-range");
        return this->_data.at(index);
    }

    /**
     * Get the cache count
     * @return Cache count
     */
    int
    count() const
    {
        return this->_count;
    }

    /**
     * Test if cache is empty
     * @return True of cache is empty otherwise false
     */
    bool
    empty() const
    {
        return (this->_count == 0);
    }

    /**
     * Test if cache is full
     * @return True if cache is full otherwise false
     */
    bool
    full() const
    {
        return (this->_count >= (int)this->_data.size());
    }

    /**
     * Reinitializes the cache sizes and clears the cache count
     */
    void
    initialize(int const size)
    {
        this->_data.resize(size, 0);
        this->reset();
    }

    /**
     * Largest cache value
     * @return Maximum value in cache
     * @throws std::exception if cache is empty
     */
    T
    max() const
    {
        if (this->empty()) throw std::runtime_error("FixedSizeCache::max() called on empty cache");
        return *std::max_element(this->_data.begin(), this->_data.begin() + this->_count);
    }

    /**
     * Smallest cache value
     * @return Minimum value in cache
     * @throws std::exception if cache is empty
     */
    T
    min() const
    {
        if (this->empty()) throw std::runtime_error("FixedSizeCache::min() called on empty cache");
        return *std::min_element(this->_data.begin(), this->_data.begin() + this->_count);
    }

    /**
     * Clears the cache but leaves the size the same
     */
    void
    reset()
    {
        this->_data.assign(this->_data.size(), 0);
        this->_count = 0;
    }

private:

    /**
     * Number of stored values in cache
     */
    int _count;

    /**
     * Cache data
     */
    std::vector<T> _data;

};

/**
 * @struct DynamicHistogramOptions
 */
struct DynamicHistogramOptions
{
    /**
     * @param[in] enabled True if histogram is enabled otherwise new data is not added histogram
     * @param[in] binCount Number of bins in histogram
     * @param[in] cacheSize Size of cache to store points not within histogram bounds
     */
    DynamicHistogramOptions(bool enabled = true, int binCount = 100, int cacheSize = 1000) :
        binCount(binCount),
        cacheSize(cacheSize),
        enabled(enabled)
    {
    }
    int binCount;
    int cacheSize;
    bool enabled;
};

/**
 * @class DynamicHistogram
 * @brief A 1D frequency histogram with uniform bins that dynamically grows as needed.
 * @details A uniform bin histogram and a cache are maintained. Any value added that is not
 *          contained within the histogram bounds is added to the cache. The histogram and
 *          cache are merged when the cache is full or when the user invokes any method that
 *          accesses the histogram. The template parameter should only be a basic numeric data
 *          type (char, short, int, float, double). The histogram is 0-indexed.
 *          Can be used as a base class where the #merge method is overridden to implement
 *          different rules for merging the histogram with the cache.
 * @warning
 *          - The histogram access methods should generally not be called before enough points
 *            have been added to the histogram. The minimum number of points before calling them
 *            is dependent on the data's statistics but should be on the order of the cache size
 *            or there is no more data to be added. If called too soon the histogram #merge function
 *            is called and may cause a poor initial estimate of the optimal bin width.
 *          - Internally all calculations are double precision, so setting the template parameter
 *            to a 64-bit integer will cast all values to a double causing rounding errors.
 */
template <class T>
class DynamicHistogram
{

public:

   /**
     * @param[in] binCount Number of bins in histogram. Must be larger than 0.
     * @param[in] cacheSize Number of cache entries to store before merging with histogram. Must be larger than 0.
     */
    DynamicHistogram(DynamicHistogramOptions const & options = DynamicHistogramOptions()) :
        _binCount(options.binCount),
        _binMax(-std::numeric_limits<double>::max()),
        _binMin(std::numeric_limits<double>::max()),
        _binWidth(0),
        _enabled(options.enabled),
        _initialized(false),
        _numMerges(0)
    {
        this->_cache.initialize(options.cacheSize);
        this->_frequencies.resize(options.binCount, 0);
        this->_temp.resize(options.binCount, 0);
    }

   /**
     */
    virtual
    ~DynamicHistogram()
    {
    }

    /**
     * Add new value to histogram.
     * @param[in] value New value to add to histogram.
     * @return True if value was added to current histogram. False if the value was added to the cache or the histogram is disabled.
     */
    bool
    add(T const value)
    {
    
        if (!this->_enabled) return false;

        double dValue = (double) value;

        /* Add new data to cache if histogram is uninitialized or not within its bounds */
        if (!(this->_initialized & this->contains_(dValue)))
        {
            /* If cache is full then merge the cache and the histogram first */
            if (this->_cache.full())
            {
                this->merge();
            }

            /* Update cache if new value is not in merged histogram */
            if (!this->contains_(dValue))
            {
                this->_cache.add(dValue);
                return false;
            }
        }

        /* Add new data to existing histogram */
        this->add_(dValue);

        return true;

    }
    
    /**
     * Bin center value for the bin located at given index.
     * @param[in] index Bin index
     * @return Bin start value.
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    double
    bin(int const index)
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::bin - Histogram could not be initialized.");
        return this->bin_(index);
    }

    /**
     * Copy of normalized cumulative distribution function.
     * @return Vector of distribution values.
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    std::vector<double>
    cdf()
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::cdf - Histogram could not be initialized.");

        std::vector<double> prob = this->pdf();

        double sum = 0;
        for (std::vector<double>::iterator it = prob.begin(); it != prob.end(); ++it)
        {
            *it += sum;
            sum = *it;
        }

        return prob;

    }

    /**
     * Test if value is contained with histogram bin bounds
     * @param[in] value Test value
     * @return True if value contained within one of its bins
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    bool
    contains(double const value)
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::contains - Histogram could not be initialized.");
        return this->contains_(value);
    }

    /**
     * Access the number of bins in the current histogram.
     * @return Number of histogram bins.
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    int
    count()
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::count - Histogram could not be initialized.");
        return this->count_();
    }

    /**
     * Histogram enabled state
     * @return True if histogram is enabled otherwise false.
     */
    bool
    enabled() const
    {
        return this->_enabled;
    }

    /**
     * The frequency value for the bin located at index.
     * @param[in] index Bin index
     * @return Frequency value.
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    double
    frequency(int const index)
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::frequency - Histogram could not be initialized.");
        return this->frequency_(index);
    }

    /**
     * Copy of frequencies.
     * @return Vector of frequency values.
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    std::vector<double>
    frequencies()
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::frequencies - Histogram could not be initialized.");
        return this->_frequencies;
    }

    /**
     * The index of the bin that contains the given value.
     * @param[in] value Source value to locate in a bin
     * @return Bin index for given value.
     * @warning If the value is out-of-range then the index is extrapolated in units of the bin width.
     *          Use #contains to ensure the value is within the histogram bin range.
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    int
    index(T const value)
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::index - Histogram could not be initialized.");
        return this->index_(value);
    }

    /**
     * Initialization state of histogram. Histogram gets initializes when cache gets full or the user invokes a getter.
     * @return True if histogram has been initialized otherwise false.
     */
    bool
    initialized() const
    {
        return this->_initialized;
    }

    /**
     * Return the maximum bin value
     * @return Maximum bin value
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    double
    max()
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::max - Histogram could not be initialized.");
        return this->_binMax;
    }

    /**
     * Return the number of histogram merges that have been performed
     * @return Number of cache merges
     */
    int
    merges()
    {
        return this->_numMerges;
    }

    /**
     * Return the minimum bin value
     * @return Minimum bin value
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    double
    min()
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::min - Histogram could not be initialized.");
        return this->_binMin;
    }

    /**
     * Estimate the k'th order statistic.
     * @details The cumulative probability of the k'th order statistic is computed as #kth/#total
     *          then the cumulative density function is inverted and linearly interpolated to
     *          estimate the k'th order statistic.
     * @return Value of k'th order statistic
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    double
    order(int kth)
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::order - Histogram could not be initialized.");

        // Estimate of the k'th order statistic
        double estimate = this->bin(0);

        // Compute order's probability and clip at [0,1]
        double numItems = this->total();
        double probability = (double) (kth) / numItems;
        probability = std::max<double>(0, probability);
        probability = std::min<double>(1, probability);

        // Look up probability in cdf
        std::vector<double> dist = this->cdf();
        int numBins = dist.size();
        for (int i = 0; i < numBins-1; i++)
        {
            if (dist.at(i+1) == 0) continue;
            if ((dist.at(i) <= probability) && (probability <= dist.at(i+1)))
            {
                if (dist.at(i+1) == dist.at(i))
                {
                    // If flat cdf avoid divide-by-zero
                    estimate = this->bin(i) + 0.5 * (this->bin(i+1) - this->bin(i));
                }
                else
                {
                    // Linearly interpolate between enclosing bins
                    double slope = (this->bin(i+1) - this->bin(i)) / (dist.at(i+1) - dist.at(i));
                    estimate = this->bin(i) + slope * (probability - dist.at(i));
                }
                return estimate;
            }
        }

        return estimate;

    }

    /**
     * Copy of normalized probability distribution function.
     * @return Vector of distribution values.
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    std::vector<double>
    pdf()
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::pdf - Histogram could not be initialized.");

        std::vector<double> prob(this->count_(),0);

        double sum = this->total();

        std::vector<double>::iterator dist = prob.begin();
        for (std::vector<double>::const_iterator freq = this->_frequencies.begin(); freq != this->_frequencies.end(); ++freq)
        {
            *dist = (*freq) / sum;
            dist++;
        }

        return prob;

    }

    /**
     * Compute the total frequency sum of the histogram.
     * @return Frequency sum.
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    double
    total()
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::total - Histogram could not be initialized.");
        double sum = std::accumulate(this->_frequencies.begin(), this->_frequencies.end(), 0.0);
        return sum;
    }

    /**
     * Bin width for the bin located at given index.
     * @param[in] index Bin index
     * @return Bin center value.
     * @throws std::runtime_exception If histogram could not be initialized.
     */
    double
    width(int const index)
    {
        if (!this->merge()) throw std::runtime_error("DynamicHistogram::width - Histogram could not be initialized.");
        return this->width_(index);
    }

    /**
     * Build a histogram of data statically. Dynamic growing is not used. However,
     * once the histogram is built it can then continue to be used dynamically.
     * @param[in] data Vector input data points to be binned.
     * @param[in] binCount Number of bins in histogram. Must be larger than 0.
     * @param[in] minValue Histogram left most minimum bin value. If null then data is pre-scanned to compute.
     * @param[in] maxValue Histogram right most maximum bin value. If null then data is pre-scanned to compute.
     * @return Output histogram with its new size equal to binCount.
     * @warning If minValue and maxValue are given, any value outside the histograms bounds will be discarded.
     */
    static
    DynamicHistogram<T>
    buildStatic(std::vector<T> & data,
                DynamicHistogramOptions const & options = DynamicHistogramOptions(),
                double * minValue = 0,
                double * maxValue = 0)
    {

        DynamicHistogram<T> histogram(options);

        // ======================================================================
        // Check if any data
        // ----------------------------------------------------------------------
        if (data.empty()) return histogram;

        // ======================================================================
        // Check if data should be scanned for min/max
        // ----------------------------------------------------------------------
        double minData = (minValue == 0 ? *std::min_element(data.begin(), data.end()) : *minValue);
        double maxData = (maxValue == 0 ? *std::max_element(data.begin(), data.end()) : *maxValue);

        // ======================================================================
        // Build exact histogram
        // ----------------------------------------------------------------------
        histogram._binMax = maxData;
        histogram._binMin = minData;
        histogram._binWidth = (histogram._binMax - histogram._binMin) / histogram._binCount;
        histogram._initialized = true;
        for (typename std::vector<T>::iterator it = data.begin(); it != data.end(); ++it)
        {
            double dValue = (double) *it;
            if (histogram.contains_(dValue)) histogram.add_(dValue);
        }

        return histogram;

    }

protected:

    /**
     * Merge the cache entries with the current histogram.
     * @return True if merge performed otherwise false
     */
    virtual
    bool
    merge()
    {

        // ======================================================================
        // Don't bother merging if no data in the cache
        // ----------------------------------------------------------------------
        if (this->_cache.empty()) return (this->_initialized ? true : false);

        // ======================================================================
        // If histogram has not been initialized then build it
        // ----------------------------------------------------------------------
        if (!this->_initialized)
        {

            // Compute new bounds
            this->_binMax = this->_cache.max();
            this->_binMin = this->_cache.min();
            this->_binWidth = (this->_binMax - this->_binMin) / (double) (this->_binCount);

            // Update with cache data
            for (int i = 0; i < this->_cache.count(); i++)
            {
                this->add_(this->_cache.at(i));
            }

            // Reset cache parameters
            this->_cache.reset();

            this->_initialized = true;

            return true;

        }

        this->_numMerges++;

        // ======================================================================
        // Determine how much to expand histogram
        // ----------------------------------------------------------------------

        // Compute dynamic range of histogram + cache data
        double cacheLeft = this->bin_(this->index_(this->_cache.min())) - 0.5 * this->width_(this->index_(this->_cache.min()));
        double cacheRight = this->bin_(this->index_(this->_cache.max())) + 0.5* this->width_(this->index_(this->_cache.max()));
        double newMin = std::min<double>(cacheLeft, this->_binMin);
        double newMax = std::max<double>(cacheRight, this->_binMax);

        // Compute how much to scale the old histogram to get the new histogram
        double oldHistogramWidth = this->_binMax - this->_binMin;
        double newHistogramWidth = newMax - newMin;
        int histogramScale = (int) ceil((newHistogramWidth / oldHistogramWidth));

        // Compute dynamic range in units of histogram scale
        int newBinCountAfterScale = (int) (histogramScale * this->_binCount);

        // Compute number of empty bins
        double newMaxAfterScale = newMin + newBinCountAfterScale * this->_binWidth;
        double emptySpace = newMaxAfterScale - newMax;
        int numEmptyBins = (int) (emptySpace / this->_binWidth);
        int numEmptyBinsToShift = (int) (numEmptyBins/2);

        // Compute final dynamic range
        double newMinFinal = newMin - numEmptyBinsToShift * this->_binWidth;
        double newMaxFinal = newMinFinal + newBinCountAfterScale * this->_binWidth;
        double newBinWidthFinal = (newMaxFinal-newMinFinal) / this->_binCount;

        // ======================================================================
        // Merge old histogram into new empty histogram
        // ----------------------------------------------------------------------
        for (int i = 0; i < this->_binCount; i++)
        {
            double oldBinCenterInNewBins = this->bin_(i);
            int newBinIndex =  (int)((oldBinCenterInNewBins - newMinFinal) / newBinWidthFinal);
            this->_temp.at(newBinIndex) += this->frequency_(i);
        }

        // ======================================================================
        // Swap new histogram with old histogram
        // ----------------------------------------------------------------------
        this->_temp.swap(this->_frequencies);
        this->_temp.assign(this->_temp.size(),0);
        this->_binMax = newMaxFinal;
        this->_binMin = newMinFinal;
        this->_binWidth = newBinWidthFinal;

        // ======================================================================
        // Add cache data to new histogram
        // ----------------------------------------------------------------------
        for (int i = 0; i < this->_cache.count(); i++)
        {
            this->add_(this->_cache.at(i));
        }
        this->_cache.reset();

        return true;

    }

private:

    /*
     * Private version of public methods that do not attempt to merge
     */
    void
    add_(double const value)
    {
        this->_frequencies.at(this->index_(value))++;
    }
    
    double
    bin_(int const index)
    {
        return this->_binMin + this->_binWidth * (index + 0.5);
    }
    
    bool
    contains_(double const value)
    {
        int idx = this->index_(value);
        bool isContained = (idx >= 0) && (idx < (int) this->_binCount);
        return isContained;
    }
    
    int
    count_()
    {
        return this->_binCount;
    }

    double
    frequency_(int const index)
    {
        return this->_frequencies.at(index);
    }

    int
    index_(double const value)
    {
        if (value == this->_binMax) return (this->_binCount - 1);
        int idx = (int) floor((value - this->_binMin) / this->_binWidth);
        return idx;
    }

    double
    width_(int const index)
    {
        return this->_binWidth;
    }

protected:

    /**
     * Number of histogram bins
     */
    int _binCount;
    
    /**
     * Histogram maximum value
     */
    double _binMax;
    
    /**
     * Histogram minimum value
     */
    double _binMin;
    
    /**
     * Histogram bin width
     */
    double _binWidth;

    /**
     * Cache where values are stored until merged into the histogram
     */
    FixedSizeCache<double> _cache;

    /**
     * Histogram initialization state
     */
    bool _enabled;

    /**
     * Histogram frequencies
     */
    std::vector<double> _frequencies;

    /**
     * Histogram initialization state
     */
    bool _initialized;

    /**
     * Number of merges performs
     */
    int _numMerges;

private:

    /**
     * Temporary container to hold new histogram frequencies during merging
     */
    std::vector<double> _temp;

};

/**
 * @struct DataPoint
 * @brief Data container structure
 */
struct DataPoint
{
    /**
     * Create default DataPoint
     */
    DataPoint() :
        active(false),
        value(0)
    {
    }

    /**
     * Create custom DataPoint
     * @param[in] value DataPoint value
     * @param[in] activate Active state of data point
     */
    DataPoint(double value,
              bool activate) :
        active(activate),
        value(value)
    {
    }

    /**
     * Active state of data point
     */
    bool active;
    
    /**
     * DataPoint value
     */
    double value;

};

/**
* @class DataVector
* @brief STL vector of #DataPoint
*/
class DataVector : public std::vector<DataPoint>
{
public:
    /**
    * @brief Set the active state of all entries to false
    */
    void
    deactivate()
    {
        for (DataVector::iterator it = this->begin(); it != this->end(); ++it)
        {
            it->active = false;
        }
    }

    /**
     * Copy subvector using mask
     */
    void
    copy(DataVector const & src, std::vector<bool> const & mask)
    {
        DataVector::copy(src, src.size(), *this, this->size(), mask, mask.size());
    }

    /**
     * Copy subvector using mask
     */
    static
    void copy(DataVector const & src,
              DataVector::size_type const & numSrc,
              DataVector & dst,
              DataVector::size_type const & numDst,
              std::vector<bool> const & mask,
              DataVector::size_type const & numMask)
    {

        // Deactivate all destination value by default
        dst.deactivate();

        // Copy elements
        DataVector::size_type j = 0;
        for (DataVector::size_type i = 0; i < numSrc; i++)
        {
            // Check if current value is masked out
            if (i < numMask)
            {
                if (!mask.at(i)) continue;
            }

            // Add or set next value
            if (j >= numDst)
            {
                dst.push_back(src.at(i));
            }
            else
            {
                dst.at(j) = src.at(i);
            }
            j++;
        }

    }
};

/**
 * @class DataFilters
 * @brief Class that tests if a data source satisfies filter requirements
 */
class DataFilters
{
    
private:

    class NumericFilter
    {
    public:
        NumericFilter(double const & minValue,
                      double const & maxValue,
                      bool const & isAccept = true) :
                        _maxValue(maxValue),
                        _minValue(minValue),
                        _isAccept(isAccept)
        {
        }
        
        /**
         * Check if token should be filtered
         */
        bool isFiltered(double const & value) const
        {
            bool isInInterval = (this->_minValue <= value) && (value <= this->_maxValue);
            return (this->_isAccept ? isInInterval : !isInInterval);
        }
        
    private:
        double _maxValue;
        double _minValue;
        bool _isAccept;
    };

    class StringFilter
    {
    public:
        StringFilter(std::string const & pattern,
                     bool const & isExact = true,
                     bool const & isSensitive = true,
                     bool const & isAccept = true) :
                        _pattern(pattern),
                        _isExact(isExact),
                        _isSensitive(isSensitive),
                        _isAccept(isAccept)
        {
            if (!isSensitive)
            {
                std::transform(this->_pattern.begin(), this->_pattern.end(), this->_pattern.begin(), ::toupper);
            }
        }
        
        /**
         * Check if token should be filtered
         */
        bool
        isFiltered(std::string const & token) const
        {

            std::string input;
            
            // Check if filtering is case sensitive
            input = token;
            if (!this->_isSensitive)
            {
                std::transform(token.begin(), token.end(), input.begin(), ::toupper);
            }
            
            // Check if filtering is exact or not
            bool isMatch = false;
            if (this->_isExact)
            {
                isMatch = input.compare(this->_pattern) == 0;
            }
            else
            {
                isMatch = input.find(this->_pattern) != std::string::npos;
            }
            
            return (this->_isAccept ? isMatch : !isMatch);
            
        }
    private:
        std::string _pattern;
        bool _isExact;
        bool _isSensitive;
        bool _isAccept;
    };

    struct NumericFilterCase
    {
        NumericFilterCase(unsigned int index, NumericFilter filter) : index(index), filter(filter)
        {
        }
        unsigned int index;
        NumericFilter filter;
    };

    struct StringFilterCase
    {
        StringFilterCase(unsigned int index, StringFilter filter) : index(index), filter(filter)
        {
        }
        unsigned int index;
        StringFilter filter;
    };
    
public:

    /*
     * Adds a numeric filter
     * @param[in] index Column index to apply filter
     * @param[in] minValue Minimum value allowed
     * @param[in] maxValue Maximum value allowed
     * @param[in] isAccept True if a match accepts value
     *                     False if a match rejects value
     */
    void
    addNumericFilter(unsigned int const & index,
                     double const & minValue,
                     double const & maxValue,
                     bool const & isAccept)
    {
        this->_numericFilters.push_back(NumericFilterCase(index, NumericFilter(minValue, maxValue, isAccept)));
    }

    /*
     * Adds a string filter
     * @param[in] index Column index to apply filter
     * @param[in] pattern Pattern to match string
     * @param[in] isExact True if string match is exact
     *                    False if string match is partial
     * @param[in] isExact True if string match is case sensitive
     *                    False if string match is case insensitive
     * @param[in] isAccept True if a match accepts value
     *                     False if a match rejects value
     */
    void
    addStringFilter(unsigned int const & index,
                    std::string const & pattern,
                    bool const & isExact,
                    bool const & isSensitive,
                    bool const & isAccept)
    {
        this->_stringFilters.push_back(StringFilterCase(index, StringFilter(pattern, isExact, isSensitive, isAccept)));
    }

    /*
     * Tests if data points passes numeric filtering
     * @param[in] data Vector of data points to test
     * @return True if a data passes filter otherwise false
     */
    bool
    isFiltered(DataVector const & data) const
    {
        unsigned int size = data.size();
        // Iterate through and check if any filter applies to this data point
        if (this->_numericFilters.size() == 0) return true;
        for (std::vector<NumericFilterCase>::const_iterator it = this->_numericFilters.begin(); it != this->_numericFilters.end(); ++it)
        {
            if (it->index >= size || !data.at(it->index).active) continue; // don't check if out-of-range or data is not active
            if (it->filter.isFiltered(data.at(it->index).value)) return true;
        }
        return false;
    }

    /*
     * Tests if strings pass string filtering
     * @param[in] data Vector of strings to test
     * @return True if a data passes filter otherwise false
     */
    bool
    isFiltered(std::vector<std::string> const & data) const
    {
        unsigned int size = data.size();
        // Iterate through and check if any filter applies to this data point
        if (this->_stringFilters.size() == 0) return true;
        for (std::vector<StringFilterCase>::const_iterator it = this->_stringFilters.begin(); it != this->_stringFilters.end(); ++it)
        {
            if (it->index >= size) continue; // don't check if out-of-range
            if (it->filter.isFiltered(data.at(it->index))) return true;
        }
        return false;
    }
    
private:

    std::vector<StringFilterCase> _stringFilters;
    std::vector<NumericFilterCase> _numericFilters;

};

/**
 * @struct StatisticsTrackerOptions
 */
struct StatisticsTrackerOptions
{
    StatisticsTrackerOptions()
    {
        this->doCov = false;
        this->doMax = false;
        this->doMean = false;
        this->doMin = false;
        this->doVar = false;
        this->histogramOptions.enabled = false;
    }
    bool doCov;
    bool doMax;
    bool doMean;
    bool doMin;
    bool doVar;
    DynamicHistogramOptions histogramOptions;
};

/**
 * @struct SamplerOptions
 */
struct SamplerOptions
{
    enum SampleMode
    {
        DEFAULT = 0,
        UNIFORM,
        RANDOM
    };

    SamplerOptions()
    {
        this->mode = SamplerOptions::DEFAULT;
        this->step = 1;
    }
    SampleMode mode;
    unsigned int step;
};

/**
 * @struct CommandLineOptions
 */
struct CommandLineOptions
{
    CommandLineOptions()
    {
        this->blankEOF = false;
        this->comment = "";
        this->delimiter = "";
        this->fileInput = "";
        this->fileOutput = "";
        this->filterColumns = std::vector<bool>();
        this->headerRow = 0;
        this->numLinesToKeep = std::numeric_limits<int>::max();
        this->numLinesToReshape = 1;
        this->numLinesToSkip = 0;
        this->removeDuplicates = false;
        this->seed = 1;
        this->showCorrelation = false;
        this->showCovariance = false;
        this->showFilteredData = false;
        this->showHistogram = false;
        this->showLeastSquaresOffset = false;
        this->showLeastSquaresSlope = false;
        this->showStatistics = false;
        this->strictParsing = false;
        this->verboseLevel = 0;
    }
    bool blankEOF;
    std::string comment;
    std::string delimiter;
    std::string fileInput;
    std::string fileOutput;
    std::vector<bool> filterColumns;
    DataFilters filterRows;
    unsigned int headerRow;
    unsigned int numLinesToKeep;
    unsigned int numLinesToReshape;
    unsigned int numLinesToSkip;
    bool removeDuplicates;
    SamplerOptions sampling;
    unsigned int seed;
    bool showCorrelation;
    bool showCovariance;
    bool showHistogram;
    bool showFilteredData;
    bool showLeastSquaresOffset;
    bool showLeastSquaresSlope;
    bool showStatistics;
    StatisticsTrackerOptions statisticsOptions;
    bool strictParsing;
    unsigned char verboseLevel;
};

/**
 * @class CommandLineParser
 */
class CommandLineParser
{

public:
    /**
     * Create command line parser
     * @param[in] argc Argument count
     * @param[in] argv Argument array
     * @throws std:runtime_error if parsing is invalid
     */
    CommandLineParser(int argc,
                      char * argv[]) :
      _showUsage(false),
      _showVersion(false)
    {
        this->parse(argc, argv);
    }

    /**
     * Print software usage
     */
    static
    void
    printUsage()
    {
        std::cout << std::endl;
        std::cout << "NAME" << std::endl;
        std::cout << std::endl;
        std::cout << "    clistats - Command line statistics tool" << std::endl;
        std::cout << std::endl;
        std::cout << "SYNOPSIS" << std::endl;
        std::cout << std::endl;
        std::cout << "    clistats [-h] [-v <level>] [-V] [-i <file>] [-o <file>] [-a <code>] [-b]" << std::endl;
        std::cout << "             [-c <character>] [-d <character>] [-fc <filter>] [-fn <filter>]" << std::endl;
        std::cout << "             [-fs <filter>] [-k <rows>] [-r] [-rs <rows>] [-s <rows>]" << std::endl;
        std::cout << "             [-se <seed>] [-sr <step>] [-st] [-su <step>] [-t <row>] [-cr]" << std::endl;
        std::cout << "             [-cv] [-fd] [-hg <bins>,[cache]] [-lo] [-ls] [-ss] < [file]" << std::endl;
        std::cout << std::endl;
        std::cout << "DESCRIPTION" << std::endl;
        std::cout << std::endl;
        std::cout << "    A command line tool to compute statistics for a stream of delimited input" << std::endl;
        std::cout << "    numbers. It takes a stream of numbers from standard input or a redirected" << std::endl;
        std::cout << "    file. To stop processing and display the statistics enter the EOF signal" << std::endl;
        std::cout << "    (CTRL-D on POSIX systems like Linux or Cygwin or CTRL-Z on Windows)." << std::endl;
        std::cout << "    Alternatively, see the --blank option to use a blank row. The default" << std::endl;
        std::cout << "    delimiter is a comma." << std::endl;
        std::cout << std::endl;
        std::cout << "    Delimited tokens that are not valid numeric values are skipped." << std::endl;
        std::cout << "    If any tokens are skipped then the resulting covariance and correlation" << std::endl;
        std::cout << "    are not exact and only a conservative overestimation." << std::endl;
        std::cout << std::endl;
        std::cout << "    Note, all statistics are computed on a rolling basis and not by using the" << std::endl;
        std::cout << "    entire dataset, so all non-order statistics are approximations." << std::endl;
        std::cout << std::endl;
        std::cout << "    Displayed statistics use headers that start with # to enable piping" << std::endl;
        std::cout << "    results to gnuplot." << std::endl;
        std::cout << std::endl;
        std::cout << "    Available statistics to display are:" << std::endl;
        std::cout << std::endl;
        std::cout << "        - Count" << std::endl;
        std::cout << "        - Minimum" << std::endl;
        std::cout << "        - Mean" << std::endl;
        std::cout << "        - Maximum" << std::endl;
        std::cout << "        - Standard deviation" << std::endl;
        std::cout << "        - Covariance" << std::endl;
        std::cout << "        - Correlation" << std::endl;
        std::cout << "        - Least Squares Slope" << std::endl;
        std::cout << "        - Least Squares Offset" << std::endl;
        std::cout << "        - Histogram" << std::endl;
        std::cout << std::endl;
        std::cout << "OPTIONS" << std::endl;
        std::cout << std::endl;
        std::cout << "    General options:" << std::endl;
        std::cout << std::endl;
        std::cout << "        -h" << std::endl;
        std::cout << "        --help" << std::endl;
        std::cout << "            Print this help message" << std::endl;
        std::cout << std::endl;
        std::cout << "        -v <level>" << std::endl;
        std::cout << "        --verbose <level>" << std::endl;
        std::cout << "            Verbose level (0 to 5). Default = 0." << std::endl;
        std::cout << std::endl;
        std::cout << "        -V" << std::endl;
        std::cout << "        --version" << std::endl;
        std::cout << "            Software version is displayed." << std::endl;
        std::cout << std::endl;
        std::cout << "    I/O options:" << std::endl;
        std::cout << std::endl;
        std::cout << "        -i <file>" << std::endl;
        std::cout << "        --input <file>" << std::endl;
        std::cout << "            Input file. Default is stdin." << std::endl;
        std::cout << std::endl;
        std::cout << "        -o <file>" << std::endl;
        std::cout << "        --output <file>" << std::endl;
        std::cout << "            Output file. Default is stdout." << std::endl;
        std::cout << std::endl;
        std::cout << "        file" << std::endl;
        std::cout << "            Input file can be redirected as standard input" << std::endl;
        std::cout << std::endl;
        std::cout << "    Parsing options:" << std::endl;
        std::cout << std::endl;
        std::cout << "        -a <code>" << std::endl;
        std::cout << "        --ascii <code>" << std::endl;
        std::cout << "            Use delimiter character based on it\'s ASCII code (e.g., 9 = tab)." << std::endl;
        std::cout << "            Can be combined with \"-d\". Multiple delimiters will be" << std::endl;
        std::cout << "            concatenated." << std::endl;
        std::cout << std::endl;
        std::cout << "        -b" << std::endl;
        std::cout << "        --blank" << std::endl;
        std::cout << "            Use blank row as the EOF signal to stop processing and display" << std::endl;
        std::cout << "            results. Default is CTRL-D for POSIX and CTRL-Z Windows." << std::endl;
        std::cout << std::endl;
        std::cout << "        -c <character>" << std::endl;
        std::cout << "        --comment <character>" << std::endl;
        std::cout << "            Input comment character. Rows starting with this character will" << std::endl;
        std::cout << "            be skipped Default is no comments are recognized." << std::endl;
        std::cout << std::endl;
        std::cout << "        -d <character>" << std::endl;
        std::cout << "        --delimiter <character>" << std::endl;
        std::cout << "            Input delimiter character. Cannot be the same as the comment" << std::endl;
        std::cout << "            character. Default = \",\". Multiple delimiters will be" << std::endl;
        std::cout << "            concatenated." << std::endl;
        std::cout << std::endl;
        std::cout << "        -fc <filter>" << std::endl;
        std::cout << "        --filterColumn" << std::endl;
        std::cout << "            Filter columns to process using comma delimited array of column" << std::endl;
        std::cout << "            indices to keep." << std::endl;
        std::cout << "            Format: \"<c1>,<c2>,...,<c3:c4>,...\"" << std::endl;
        std::cout << "                (e.g., \"1,2,3,5\" will discard column 4). Range of" << std::endl;
        std::cout << "                columns can be indicated using \":\" (e.g., 1:3,5)." << std::endl;
        std::cout << std::endl;
        std::cout << "        -fn <filter>" << std::endl;
        std::cout << "        --filterNumeric <filter>" << std::endl;
        std::cout << "            Filter rows to process only those with a column matching a numeric" << std::endl;
        std::cout << "            criteria. Multiple filters can be used and will be processed in the" << std::endl;
        std::cout << "            order the -fn flags are provided." << std::endl;
        std::cout << "            Format: \"#,min,max,[a|r]\"" << std::endl;
        std::cout << "                # = Number corresponding to column index to filter" << std::endl;
        std::cout << "                min = Minimum numeric value (-inf for no minimum)" << std::endl;
        std::cout << "                max = Maximum numeric value (inf for no maximum)" << std::endl;
        std::cout << "                a|r = Optional letter to denote accept/reject filtering" << std::endl;
        std::cout << "                      \"a\" denotes row is accepted if criteria is met" << std::endl;
        std::cout << "                      \"r\" denotes row is rejected if criteria is met" << std::endl;
        std::cout << std::endl;
        std::cout << "        -fs <filter>" << std::endl;
        std::cout << "        --filterString <filter>" << std::endl;
        std::cout << "            Filter rows to process only those with a column matching a string" << std::endl;
        std::cout << "            match. Multiple filters can be used and will be processed in the" << std::endl;
        std::cout << "            order the -fs flags are provided." << std::endl;
        std::cout << "            Format: \"#,p,[e|p],[s|i],[a|r]\"" << std::endl;
        std::cout << "                # = Number corresponding to column index to filter" << std::endl;
        std::cout << "                p = String of pattern to match" << std::endl;
        std::cout << "                e|p = Optional letter to denote exact or partial match" << std::endl;
        std::cout << "                s|i = Optional letter to denote case sensitive/insensitive" << std::endl;
        std::cout << "                a|r = Optional letter to denote accept/reject filtering" << std::endl;
        std::cout << "                      \"a\" denotes row is accepted if criteria is met" << std::endl;
        std::cout << "                      \"r\" denotes row is rejected if criteria is met" << std::endl;
        std::cout << std::endl;
        std::cout << "        -k <rows>" << std::endl;
        std::cout << "        --keep <rows>" << std::endl;
        std::cout << "            Number of rows to keep after processing begins." << std::endl;
        std::cout << "            Default is to keep all." << std::endl;
        std::cout << std::endl;
        std::cout << "        -r" << std::endl;
        std::cout << "        --remove" << std::endl;
        std::cout << "            Remove duplicate sequential delimiters and treat them as one" << std::endl;
        std::cout << "            delimiter. Default is to treat them separately." << std::endl;
        std::cout << std::endl;
        std::cout << "        -rs <rows>" << std::endl;
        std::cout << "        --reshape <rows>" << std::endl;
        std::cout << "            Reshape every N rows in a single column. Default = 1." << std::endl;
        std::cout << std::endl;
        std::cout << "        -s <rows>" << std::endl;
        std::cout << "        --skip <rows>" << std::endl;
        std::cout << "            Number of rows to skip before processing begins (e.g., header row)" << std::endl;
        std::cout << "            Default = 0 (process all rows)." << std::endl;
        std::cout << std::endl;
        std::cout << "        -se <seed>" << std::endl;
        std::cout << "        --seed <seed>" << std::endl;
        std::cout << "            Seed for random number generator when using \"-sr\" option." << std::endl;
        std::cout << "            <seed> can take values in the range [1," << RAND_MAX << "]." << std::endl;
        std::cout << "            If set to 0 then the current time is used as the seed." << std::endl;
        std::cout << "            Default = 1." << std::endl;
        std::cout << std::endl;
        std::cout << "        -sr <step>" << std::endl;
        std::cout << "        --samplerandom <step>" << std::endl;
        std::cout << "            Sample the rows to be processed randomly in intervals of given" << std::endl;
        std::cout << "            step size. One row is randomly chosen every \"step\" rows." << std::endl;
        std::cout << "            Use \"-se\" to change the random number generator seed." << std::endl;
        std::cout << "            Mutually exclusive with option \"-su\"." << std::endl;
        std::cout << std::endl;
        std::cout << "        -st" << std::endl;
        std::cout << "        --strict" << std::endl;
        std::cout << "            Strictly enforce parsing to only accept rows with the same number" << std::endl;
        std::cout << "            of tokens as the first successfully parsed row. Default is to" << std::endl;
        std::cout << "            allow missing tokens." << std::endl;
        std::cout << std::endl;
        std::cout << "        -su <step>" << std::endl;
        std::cout << "        --sampleuniform <step>" << std::endl;
        std::cout << "            Sample the rows to be processed uniformly in intervals of given" << std::endl;
        std::cout << "            step size. Mutually exclusive with option \"-sr\"." << std::endl;
        std::cout << std::endl;
        std::cout << "        -t <row>" << std::endl;
        std::cout << "        --titles <row>" << std::endl;
        std::cout << "            Row containing column header titles. If absent then index is used" << std::endl;
        std::cout << "            for each column title." << std::endl;
        std::cout << std::endl;
        std::cout << "    Display options:" << std::endl;
        std::cout << std::endl;
        std::cout << "        -cr" << std::endl;
        std::cout << "        --correlation" << std::endl;
        std::cout << "            Show correlation of data." << std::endl;
        std::cout << std::endl;
        std::cout << "        -cv" << std::endl;
        std::cout << "        --covariance" << std::endl;
        std::cout << "            Show covariance of data." << std::endl;
        std::cout << std::endl;
        std::cout << "        -fd" << std::endl;
        std::cout << "        --filtered" << std::endl;
        std::cout << "            Show filtered data to be processed." << std::endl;
        std::cout << std::endl;
        std::cout << "        -hg <bins>,[cache]" << std::endl;
        std::cout << "        --histogram <bins>,[cache]" << std::endl;
        std::cout << "            Compute the histogram. Number of bins is required (Default = 100)." << std::endl;
        std::cout << "            Cache size is optional and must be greater than 100" << std::endl;
        std::cout << "            (Default = 1000)." << std::endl;
        std::cout << std::endl;
        std::cout << "        -lo" << std::endl;
        std::cout << "        --lsqoffset" << std::endl;
        std::cout << "            Show least squares linear offset." << std::endl;
        std::cout << std::endl;
        std::cout << "        -ls" << std::endl;
        std::cout << "        --lsqslope" << std::endl;
        std::cout << "            Show least squares linear slope." << std::endl;
        std::cout << std::endl;
        std::cout << "        -ss" << std::endl;
        std::cout << "        --statistics" << std::endl;
        std::cout << "            Show summary statistics. Enabled by default unless another display" << std::endl;
        std::cout << "            option is enabled." << std::endl;
        std::cout << std::endl;
        std::cout << "EXIT STATUS" << std::endl;
        std::cout << std::endl;
        std::cout << "    0 = Success" << std::endl;
        std::cout << "    1 = Argument parsing failed" << std::endl;
        std::cout << "    2 = Computing statistics failed" << std::endl;
        std::cout << "    3 = Displaying statistics failed" << std::endl;
        std::cout << std::endl;
        std::cout << "BUGS" << std::endl;
        std::cout << std::endl;
        std::cout << "    Report bugs to " << ApplicationProperties::AUTHOR << std::endl;
    }

    /**
     * Print software version
     */
    static
    void
    printVersion()
    {
        std::cout << ApplicationProperties::version() << std::endl;
    }

    /**
     * Application usage flag
     * @return True if usage was requested otherwise false
     */
    bool
    showUsage() const
    {
        return this->_showUsage;
    }

    /**
     * Application version flag
     * @return True if version was requested otherwise false
     */
    bool
    showVersion() const
    {
        return this->_showVersion;
    }
    
private:

    /**
     * Parse command line arguments
     * @param[in] argc Argument count
     * @param[in] argv Argument array
     * @throws std::runtime_error if flag argument is unrecognized or parsing its value failed
     */
    void
    parse(int argc,
          char * argv[])
    {

        // Iterate through each argument until all are parsed
        // If argument is recognized as a valid flag then it is parsed
        //      If parsing fails on flag value an exception is thrown
        // If argument is not recognized as a valid flag an exception is thrown
        std::list<std::string> arguments(argv, argv + argc);
        arguments.pop_front(); // Remove executable name
        while (!arguments.empty())
        {
            std::list<std::string>::iterator argument = arguments.begin();
            std::string flag(*argument);
            if (flag == "-h" || flag == "--help")
            {
                this->_showUsage = true;
                return;
            }
            else if (flag == "-V" || flag == "--version")
            {
                this->_showVersion = true;
                return;
            }
            else if (flag == "-v" || flag == "--verbose")
            {
                // Parse verbosity flag value for a single integer within the range of log levels
                unsigned int verbose = 0;
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                bool isInt = StringParser::toValue<unsigned int>(value, verbose);
                if (!(isInt && verbose >= Logger::Level::FATAL && verbose <= Logger::Level::DETAIL))
                {
                    throw std::runtime_error("Invalid verbose level");
                }
                this->options.verboseLevel = (unsigned char)verbose;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-i" || flag == "--input")
            {
                // Parse input file flag value for a single string
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                std::ifstream f(value.c_str());
                if (!f.good())
                {
                    throw std::runtime_error("Invalid input file");
                }
                this->options.fileInput = value;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-o" || flag == "--output")
            {
                // Parse output file flag value for a single string
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                std::ofstream f(value.c_str());
                if (!f.good())
                {
                    throw std::runtime_error("Invalid output file");
                }
                this->options.fileOutput = value;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-a" || flag == "--ascii")
            {
                // Parse delimiter flag value for a ASCII code character
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                // Covert code to character
                int code = 0;
                if (!(StringParser::toValue<int>(value,code) && (code >= 0) && (code <= 127)))
                {
                    throw std::runtime_error("Invalid delimiter ASCII code");
                }
                char character = (char) code;
                this->options.delimiter += character;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-b" || flag == "--blank")
            {
                this->options.blankEOF = true;
                arguments.pop_front();
            }
            else if (flag == "-c" || flag == "--comment")
            {
                // Parse comment flag value for a single character
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                if (value.size() != 1)
                {
                    throw std::runtime_error("Invalid comment character = " + value);
                }
                this->options.comment = value;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-d" || flag == "--delimiter")
            {
                // Parse delimiter flag value for a single character
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                this->options.delimiter += value;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-s" || flag == "--skip")
            {
                // Parse skip flag value for a single non-negative integer
                int skip = 0;
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                bool isInt = StringParser::toValue<int>(value, skip);
                if (!(isInt && skip >= 0))
                {
                    throw std::runtime_error("Invalid number of lines to skip");
                }
                this->options.numLinesToSkip = skip;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-k" || flag == "--keep")
            {
                // Parse keep flag value for a single non-negative integer
                int keep = 0;
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                bool isInt = StringParser::toValue<int>(value, keep);
                if (!(isInt && keep >= 0))
                {
                    throw std::runtime_error("Invalid number of lines to keep");
                }
                this->options.numLinesToKeep = keep;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-r" || flag == "--remove")
            {
                this->options.removeDuplicates = true;
                arguments.pop_front();
            }
            else if (flag == "-t" || flag == "--titles")
            {
                // Parse titles flag value for a single positive integer
                unsigned int row = 0;
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                bool isInt = StringParser::toValue<unsigned int>(value, row);
                if (!(isInt && row > 0))
                {
                    throw std::runtime_error("Invalid header title row");
                }
                this->options.headerRow = row;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-fc" || flag == "--filterColumn")
            {
                // Parse filter flag value for a comma delimited vector of positive integers
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                
                // Parse filters into delimited tokens
                StringSplitter parser(value, ",");
                std::vector<std::string> values = parser.tokens();
                
                std::vector<int> filters;
                for (unsigned int j = 0; j < values.size(); j++)
                {
                    // Parse string to integers and check for ranges (e.g., 3:5)
                    std::vector<int> range;
                    bool isRangeValid = StringSplitter::toIntegers(values.at(j), range);
                    if (!isRangeValid)
                    {
                        throw std::runtime_error("Invalid filter column value = " + values.at(j));
                    }
                    
                    // Check each range entry is positive
                    for (std::vector<int>::iterator it = range.begin(); it != range.end(); ++it)
                    {
                        if (*it <= 0)
                        {
                            throw std::runtime_error("Invalid filter column value = " + values.at(j));
                        }
                        filters.push_back(*it-1);
                    }                    
                    
                }
                
                // Populate filter mask
                int maxFilter = *std::max_element(filters.begin(), filters.end());
                this->options.filterColumns.resize(maxFilter+1, false);
                for (unsigned int j = 0; j < filters.size(); j++)
                {
                    this->options.filterColumns.at(filters.at(j)) = true;
                }
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-fn" || flag == "--filterNumeric")
            {
                // Parse filter flag value for a comma delimited vector of strings
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                
                // Parse filters into delimited tokens
                StringSplitter parser(value, ",");
                std::vector<std::string> values = parser.tokens();
                
                if (values.size() == 0)
                {
                    throw std::runtime_error("Invalid numeric filter value = " + value);
                }

                // Numeric filter
                if (values.size() < 3 || values.size() > 4)
                {
                    throw std::runtime_error("Invalid numeric filter format = " + value);
                }
                
                // Check filter index is a positive integer
                int column = 0;
                bool isInt = StringParser::toValue<int>(values.at(0), column);
                if (!(isInt && column > 0))
                {
                    throw std::runtime_error("Invalid numeric filter column value = " + values.at(0));
                }
                
                // Check minimum value is a positive integer
                double minValue = 0;
                if (values.at(1).compare("-inf")==0)
                {
                    minValue = -std::numeric_limits<double>::max();
                }
                else if (values.at(1).compare("inf")==0)
                {
                    minValue = std::numeric_limits<double>::max();
                }
                else
                {
                    bool isDouble = StringParser::toValue<double>(values.at(1), minValue);
                    if (!isDouble)
                    {
                        throw std::runtime_error("Invalid numeric filter minimum value = " + values.at(1));
                    }
                }
                
                // Check maximum value is a positive integer
                double maxValue = 0;
                if (values.at(2).compare("-inf")==0)
                {
                    maxValue = -std::numeric_limits<double>::max();
                }
                else if (values.at(2).compare("inf")==0)
                {
                    maxValue = std::numeric_limits<double>::max();
                }
                else
                {
                    bool isDouble = StringParser::toValue<double>(values.at(2), maxValue);
                    if (!isDouble)
                    {
                        throw std::runtime_error("Invalid numeric filter maximum value = " + values.at(2));
                    }
                }
                
                // Check for optional filter matching
                bool isAccept = true;
                if (values.size() >= 4)
                {
                    std::string accept = values.at(3);
                    if ((accept.compare("a") != 0) && (accept.compare("r") != 0))
                    {
                        throw std::runtime_error("Invalid string filter accept option = " + values.at(3));
                    }
                    isAccept = (accept.compare("a") == 0);
                }

                // Build filter
                this->options.filterRows.addNumericFilter(column-1, minValue, maxValue, isAccept);

                arguments.pop_front();
                arguments.pop_front();
                
            }
            else if (flag == "-fs" || flag == "--filterString")
            {
                // Parse filter flag value for a comma delimited vector of strings
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                
                // Parse filters into delimited tokens
                StringSplitter parser(value, ",");
                std::vector<std::string> values = parser.tokens();
                
                // Check string filter
                if (values.size() < 2 || values.size() > 5)
                {
                    throw std::runtime_error("Invalid string filter format = " + value);
                }
                
                // Check filter index is a positive integer
                int column = 0;
                bool isInt = StringParser::toValue<int>(values.at(0), column);
                if (!(isInt && column > 0))
                {
                    throw std::runtime_error("Invalid string filter column value = " + values.at(0));
                }
                
                // Get filter pattern
                std::string pattern = values.at(1);
                
                // Check filter matching option
                bool isExact = true;
                if (values.size() >= 3)
                {
                    std::string match = values.at(2);
                    if ((match.compare("e") != 0) && (match.compare("p") != 0))
                    {
                        throw std::runtime_error("Invalid string filter matching option = " + values.at(2));
                    }
                    isExact = (match.compare("e") == 0);
                }

                // Check filter case sensitivity option
                bool isSensitive = true;
                if (values.size() >= 4)
                {
                    std::string sensitivity = values.at(3);
                    if ((sensitivity.compare("s") != 0) && (sensitivity.compare("i") != 0))
                    {
                        throw std::runtime_error("Invalid string filter case sensitivity option = " + values.at(3));
                    }
                    isSensitive = (sensitivity.compare("s") == 0);
                }
                
                // Check filter matching option
                bool isAccept = true;
                if (values.size() >= 5)
                {
                    std::string accept = values.at(4);
                    if ((accept.compare("a") != 0) && (accept.compare("r") != 0))
                    {
                        throw std::runtime_error("Invalid string filter accept option = " + values.at(4));
                    }
                    isAccept = (accept.compare("a") == 0);
                }

                // Build filter
                this->options.filterRows.addStringFilter(column-1, pattern, isExact, isSensitive, isAccept);

                arguments.pop_front();
                arguments.pop_front();
                
            }
            else if (flag == "-rs" || flag == "--reshape")
            {
                // Parse skip flag value for a single positive integer
                int reshape = 0;
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                bool isInt = StringParser::toValue<int>(value, reshape);
                if (!(isInt && reshape > 0))
                {
                    throw std::runtime_error("Invalid number of lines to reshape");
                }
                this->options.numLinesToReshape = reshape;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-se" || flag == "--seed")
            {
                // Parse skip flag value for a single positive integer
                int seed = 0;
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                bool isInt = StringParser::toValue<int>(value, seed);
                if (!(isInt && seed >= 0))
                {
                    throw std::runtime_error("Invalid seed value");
                }
                this->options.seed = seed;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-sr" || flag == "--samplerandom")
            {
                // Parse step flag value for a single non-negative integer
                int step = 0;
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                bool isInt = StringParser::toValue<int>(value, step);
                if (!(isInt && step >= 1))
                {
                    throw std::runtime_error("Invalid random sampling step size. Must be >= 1.");
                }
                if (this->options.sampling.mode != SamplerOptions::DEFAULT)
                {
                    throw std::runtime_error("Cannot use multiple sampling options.");
                }
                this->options.sampling.mode = SamplerOptions::RANDOM;
                this->options.sampling.step = step;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-st" || flag == "--strict")
            {
                this->options.strictParsing = true;
                arguments.pop_front();
            }
            else if (flag == "-su" || flag == "--sampleuniform")
            {
                // Parse step flag value for a single non-negative integer
                int step = 0;
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }
                bool isInt = StringParser::toValue<int>(value, step);
                if (!(isInt && step >= 1))
                {
                    throw std::runtime_error("Invalid uniform sampling step size. Must be >= 1.");
                }
                if (this->options.sampling.mode != SamplerOptions::DEFAULT)
                {
                    throw std::runtime_error("Cannot use multiple sampling options.");
                }
                this->options.sampling.mode = SamplerOptions::UNIFORM;
                this->options.sampling.step = step;
                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-cv" || flag == "--covariance")
            {
                this->options.showCovariance = true;
                arguments.pop_front();
            }
            else if (flag == "-cr" || flag == "--correlation")
            {
                this->options.showCorrelation = true;
                arguments.pop_front();
            }
            else if (flag == "-fd" || flag == "--filtered")
            {
                this->options.showFilteredData = true;
                arguments.pop_front();
            }
            else if (flag == "-hg" || flag == "--histogram")
            {
                // Parse histogram flag value for a single positive integer
                this->options.showHistogram = true;
                int binCount = 0;
                std::string value = "";
                if ((++argument) != arguments.end())
                {
                    value = *argument;
                }

                // Parse filters into delimited tokens
                StringSplitter parser(value, ",");
                std::vector<std::string> values = parser.tokens();

                if (values.size() == 0 || (values.size() == 1 && values.at(0).empty()))
                {
                    throw std::runtime_error("Missing histogram options");
                }

                if (values.size() > 2)
                {
                    throw std::runtime_error("Invalid histogram options = " + value);
                }

                // Parse bin count
                bool isInt = StringParser::toValue<int>(values.at(0), binCount);
                if (!(isInt && binCount > 0))
                {
                    throw std::runtime_error("Invalid number of bins = " + values.at(0));
                }
                this->options.statisticsOptions.histogramOptions.binCount = binCount;

                // Parse cache size
                if (values.size() == 2)
                {
                    int cacheSize = 0;
                    bool isInt = StringParser::toValue<int>(values.at(1), cacheSize);
                    if (!(isInt && cacheSize >= 100))
                    {
                        throw std::runtime_error("Invalid histogram cache size = " + values.at(1));
                    }
                    this->options.statisticsOptions.histogramOptions.cacheSize = cacheSize;
                }

                arguments.pop_front();
                arguments.pop_front();
            }
            else if (flag == "-lo" || flag == "--lsqoffset")
            {
                this->options.showLeastSquaresOffset = true;
                arguments.pop_front();
            }
            else if (flag == "-ls" || flag == "--lsqslope")
            {
                this->options.showLeastSquaresSlope = true;
                arguments.pop_front();
            }
            else if (flag == "-ss" || flag == "--statistics")
            {
                this->options.showStatistics = true;
                arguments.pop_front();
            }
            else
            {
                throw std::runtime_error("Unrecognized flag (" + flag + "). Run \"clistats --help\" for usage.");
            }

        }

        // Check if delimiter has been provided
        if (this->options.delimiter.empty())
        {
            this->options.delimiter = ",";
        }

        // Check for conflicting options
        if (this->options.delimiter == this->options.comment)
        {
            throw std::runtime_error("Delimiter and comment characters cannot be the same.");
        }

        // Check if default basic statistics should be shown
        if (!this->options.showStatistics)
        {
            this->options.showStatistics = !(this->options.showCovariance ||
                                             this->options.showCorrelation ||
                                             this->options.showFilteredData ||
                                             this->options.showHistogram ||
                                             this->options.showLeastSquaresOffset ||
                                             this->options.showLeastSquaresSlope);
        }

        // Determine what statistics should be computed
        this->options.statisticsOptions.doCov = (this->options.showCovariance ||
                                                 this->options.showCorrelation ||
                                                 this->options.showLeastSquaresOffset ||
                                                 this->options.showLeastSquaresSlope);
        this->options.statisticsOptions.doMax = this->options.showStatistics;
        this->options.statisticsOptions.doMean = this->options.showStatistics || this->options.statisticsOptions.doCov;
        this->options.statisticsOptions.doMin = this->options.showStatistics;
        this->options.statisticsOptions.doVar = this->options.showStatistics;
        this->options.statisticsOptions.histogramOptions.enabled = this->options.showHistogram;

    }

public:

    /**
     * Options
     */
    CommandLineOptions options;

private:

    /**
     * Show application usage flag
     */
    bool _showUsage;

    /**
     * Show application version flag
     */
    bool _showVersion;

};

/**
 * @class StatisticsTracker
 * @brief Tracks the statistics of a scalar variable.
 * @details The order determines which statistics will be tracked.
 *            - 0 = 0th order statistic: minimum, maximum
 *            - 1 = 1st order statistic: mean
 *            - 2 = 2nd order statistic: variance
 *            - >3 = Not supported.
 */
class StatisticsTracker
{

public:

    /**
     * Custom constructor
     * @param[in] name Statistic name
     * @param[in] count Initial number of data points
     * @param[in] minimum Initial data minimum
     * @param[in] mean Initial data mean
     * @param[in] maximum Initial data maximum
     * @param[in] variance Initial data variance
     * @param[in] options Statistics options
     */
    StatisticsTracker(std::string name = "",
                      long count = 0,
                      double minimum = 0,
                      double mean = 0,
                      double maximum = 0,
                      double variance = 0,
                      StatisticsTrackerOptions const & options = StatisticsTrackerOptions()) :
        _count(count),
        _histogram(DynamicHistogram<double>(options.histogramOptions)),
        _maximum(maximum),
        _mean(mean),
        _minimum(minimum),
        _name(name),
        _options(options),
        _variance(variance)
    {
    }

    /**
     * Get data count
     * @return Number of data points
     */
    long
    getCount(void) const
    {
        return this->_count;
    }

    /**
     * Get data histogram
     * @return Histogram
     */
    DynamicHistogram<double>
    getHistogram(void) const
    {
        return this->_histogram;
    }
    
    /**
     * Get data minimum
     * @return Data minimum
     */
    double
    getMinimum(void) const
    {
        return this->_minimum;
    }

    /**
     * Get data mean
     * @return Data mean
     */
    double
    getMean(void) const
    {
        return this->_mean;
    }

    /**
     * Get data maximum
     * @return Data maximum
     */
    double
    getMaximum(void) const
    {
        return this->_maximum;
    }

    /**
     * Get data name
     * @return Data name
     */
    std::string
    getName(void) const
    {
        return this->_name;
    }

    /**
     * Get data variance
     * @return Data variance
     */
    double
    getVariance(void) const
    {
        return this->_variance;
    }

    /**
     * Set data name
     * @param[in] name Data name
     */
    void
    setName(std::string const & name)
    {
        this->_name = name;
    }

    /**
     * Update statistics with new data value
     * @param[in] value Data value
     */
    void
    update(double const value)
    {

        // update counter
        this->_count++;

        // initialize if first value
        if (this->_count == 1)
        {
            if (this->_options.doMin)
            {
                this->_minimum = value;
            }
            if (this->_options.doMax)
            {
                this->_maximum = value;
            }
            if (this->_options.doMean)
            {
                this->_mean = value;
            }
            if (this->_options.doVar)
            {
                this->_variance = 0;
            }
            if (this->_options.histogramOptions.enabled)
            {
                this->_histogram.add(value);
            }
            return;
        }

        // update min
        if (this->_options.doMin)
        {
            if (value < this->_minimum)
            {
                this->_minimum = value;
            }
        }

        // update max
        if (this->_options.doMax)
        {
            if (value > this->_maximum)
            {
                this->_maximum = value;
            }
        }

        if (this->_options.doMean || this->_options.doVar)
        {
            double delta = value - this->_mean;

            if (this->_options.doMean)
            {
                // update mean
                this->_mean += delta / (double) this->_count;
            }

            if (this->_options.doVar)
            {
                // update variance
                double scale = (this->_count - 1.0) / (double) this->_count;
                this->_variance *= scale;
                this->_variance += delta * delta * scale / (double) this->_count;
                
            }

        }

        if (this->_options.histogramOptions.enabled)
        {
            this->_histogram.add(value);
        }

    }

protected:

    /**
     * Statistic count
     */
    long _count;

    /**
     * Statistic histogram
     */
    DynamicHistogram<double> _histogram;

    /**
     * Statistic maximum
     */
    double _maximum;

    /**
     * Statistic mean
     */
    double _mean;

    /**
     * Statistic minimum
     */
    double _minimum;

    /**
     * Statistic name
     */
    std::string _name;
    
    /**
     * Statistics options
     */
    StatisticsTrackerOptions _options;

    /**
     * Statistic variance
     */
    double _variance;

};

/**
 * @class MultivariateTracker
 * @brief Tracks the statistics of a multivariate variable.
 * @details   The order determines which statistics will be tracked.
 *              - 0 = 0th order statistic: minimum, maximum
 *              - 1 = 1st order statistic: mean
 *              - 2 = 2nd order statistic: variance
 *              - >3 = Not supported.
 */
class MultivariateTracker
{

public:

    /**
     * Custom constructor
     * @param[in] size Dimensional size of statistics to track
     * @param[in] options Statistics options
     */
    MultivariateTracker(long size = 0,
                        StatisticsTrackerOptions const & options = StatisticsTrackerOptions()) :
        _covariance(std::vector<double>(size*size,0)),
        _options(options)
    {
        this->initialize(size, options);
    }

    /**
     * Get data count
     * @param[in] index Index of multi-variable
     * @return Number of data points
     */
    long
    getCount(long const index) const
    {
        return this->_statistics[index].getCount();
    }

    /**
     * Get data covariance
     * @param[in] i Index of i'th variable
     * @param[in] j Index of j'th variable
     * @return Data covariance
     * @throws std::runtime_error if invalid indices
     */
    double
    getCovariance(long const i,
                  long const j) const
    {
        const int numPoints = this->_statistics.size();
        if ((i < 0) || (i >= numPoints) || (j < 0) || (j >= numPoints)) throw std::runtime_error("Invalid covariance indices");
        int k = i * numPoints + j;
        return this->_covariance.at(k);
    }

    /**
     * Get data correlation
     * @param[in] i Index of i'th variable
     * @param[in] j Index of j'th variable
     * @return Data correlation
     * @throws std::runtime_error if invalid indices
     */
    double
    getCorrelation(long const i,
                   long const j) const
    {
        double varI = this->getCovariance(i,i);
        double varJ = this->getCovariance(j,j);
        return this->getCovariance(i,j) / (sqrt(varI) * sqrt(varJ));
    }

    /**
     * Get tracker dimension
     * @return Tracker dimension
     */
    long
    getDimension() const
    {
        return (long) this->_statistics.size();
    }

    /**
     * Get data histogram
     * @param[in] index Index of multi-variable
     * @return Data histogram
     */
    DynamicHistogram<double>
    getHistogram(long const index) const
    {
        return this->_statistics[index].getHistogram();
    }

    /**
     * Get data minimum
     * @param[in] index Index of multi-variable
     * @return Data minimum
     */
    double
    getMinimum(long const index) const
    {
        return this->_statistics[index].getMinimum();
    }

    /**
     * Get data mean
     * @param[in] index Index of multi-variable
     * @return Data mean
     */
    double
    getMean(long const index) const
    {
        return this->_statistics[index].getMean();
    }

    /**
     * Get data maximum
     * @param[in] index Index of multi-variable
     * @return Data maximum
     */
    double
    getMaximum(long const index) const
    {
        return this->_statistics[index].getMaximum();
    }

    /**
     * Get data name
     * @param[in] index Index of multi-variable
     * @return Data name
     */
    std::string
    getName(long const index) const
    {
        return this->_statistics[index].getName();
    }

    /**
     * Get least squares offset between i'th (independent) and j'th (independent) variables
     * @param[in] i Index of i'th variable
     * @param[in] j Index of j'th variable
     * @return Least squares offset
     * @throws std::runtime_error if invalid indices
     */
    double
    getLeastSquaresOffset(long const i,
                          long const j) const
    {
        return this->getMean(j) -  this->getLeastSquaresSlope(i,j) * this->getMean(i);
    }

    /**
     * Get least squares slope between i'th (independent) and j'th (independent) variables
     * @param[in] i Index of i'th variable
     * @param[in] j Index of j'th variable
     * @return Least squares slope
     * @throws std::runtime_error if invalid indices
     */
    double
    getLeastSquaresSlope(long const i,
                         long const j) const
    {
        return this->getCovariance(i,j) / this->getCovariance(i,i);
    }

    /**
     * Get data variance
     * @param[in] index Index of multi-variable
     * @return Data variation
     */
    double
    getVariance(long const index) const
    {
        return this->_statistics[index].getVariance();
    }

    /**
     * Set data name
     * @param[in] index Index of multi-variable
     * @param[in] name Data name
     */
    void
    setName(long const index,
            std::string name)
    {
        this->_statistics[index].setName(name);
    }

    /**
     * Update statistics with new data points
     * @param[in] data New data point
     * @return True if statisticc updated otherwise false
     */
    bool
    update(DataVector const & data)
    {

        // Check data dimensions matches tracker dimensions
        if (this->_statistics.size() != 0 && data.size() != this->_statistics.size())
        {
            return false;
        }

        // Try to initialize the tracker
        if (this->_statistics.size() == 0)
        {
            if (data.size() == 0) return false;
            this->initialize(data.size(), this->_options);
        }

        if (this->_options.doCov)
        {

            // Update covariance
            const int numPoints = data.size();
            for (int i = 0; i < numPoints; i++)
            {
            
                if (!data[i].active) continue;

                int countI = this->_statistics[i].getCount();
                double meanI = this->_statistics[i].getMean();
                double deltaI = (data[i].value - meanI);

                for (int j = 0; j < numPoints; j++)
                {

                    if (!data[j].active) continue;
                
                    int countJ = this->_statistics[j].getCount();
                    double meanJ = this->_statistics[j].getMean();
                    double deltaJ = (data[j].value - meanJ);

                    // Conservative estimate the number of data points
                    double N = std::max<double>(countI+1, countJ+1);

                    int k = i * numPoints + j;
                    this->_covariance[k] *= (N-1);
                    this->_covariance[k] += ((N-1)/N) * (deltaI) * (deltaJ);
                    this->_covariance[k] /= N;

                }
            }

        }

        // Update statistics with each data point
        DataVector::const_iterator datum = data.begin();
        std::vector<StatisticsTracker>::iterator tracker = this->_statistics.begin();
        for ( ; tracker != this->_statistics.end(); ++tracker, ++datum)
        {
            if (datum->active) tracker->update(datum->value);
        }

        return true;
        
    }

private:

    /**
     * Initialize the statistics array
     * @param[in] size Size of statistics array
     * @param[in] options Statistics options
     */
    void
    initialize(long size,
               StatisticsTrackerOptions const & options)
    {
        // Initialize scalar statistics trackers
        for (long i = 0 ; i < size; ++i)
        {
            this->_statistics.push_back(StatisticsTracker("", 0, 0, 0, 0, 0, options));
        }
        
        // Initialize multivariate statistics
        this->_covariance.resize(size*size,0);
    }

protected:

    /**
     * Multivariate statistics
     */
    std::vector<double> _covariance;
    
    /**
     * Statistics options
     */
    StatisticsTrackerOptions _options;

    /**
     * Statistic trackers
     */
    std::vector<StatisticsTracker> _statistics;
    
};

/**
 * @class StatisticsWriter
 * @brief Static class to write various statistics to a stream.
 */
class StatisticsWriter
{
public:

    /**
     * @class StatisticsWriter
     * @brief Static class to write various statistics to a stream.
     * @param[in] comment Character used to comment header rows. Default "#".
     */
    StatisticsWriter(std::string comment = "#") :
        _comment(comment)
    {
    }

    /**
     * Write data vector to stream
     * @param[out] stream Destination stream
     * @param[in] data Source data vector to write out
     * @param[in] delimiter Data output delimiter
     */
    void
    writeData(std::ostream & stream,
              DataVector const & data,
              std::string const delimiter = " ")
    {

        for (DataVector::const_iterator it = data.begin(); it != data.end();)
        {
            // Skip inactive data
            if (!it->active)
            {
                ++it;
                continue;
            }

            // Write out data
            stream << it->value;
            ++it;

            // Write out delimiter if not last entry
            if (it != data.end())
            {
                stream << delimiter;
            }
        }
        stream << std::endl;

    }

    /**
     * Write tracker histograms to a stream
     * @param[out] stream Destination stream
     * @param[in] tracker Source tracker
     */
    void
    writeHistograms(std::ostream & stream,
                    MultivariateTracker & tracker)
    {
        
        int numDim = tracker.getDimension();

        // Print no data header
        if (numDim == 0)
        {
            int totalWidth = 20;
            std::string header(totalWidth,'='); 
            stream << this->_comment << header << std::endl;
            stream << this->_comment << std::setfill (' ') << std::setw (totalWidth/2) << "Histograms" << std::endl;
            stream << this->_comment << header << std::endl;
            stream << this->_comment << "No valid data" << std::endl;
            return;
        }

        for (long i = 0; i < numDim; i++)
        {
            const std::string name = tracker.getName(i);
            DynamicHistogram<double> histogram = tracker.getHistogram(i);
            this->writeHistogram<double>(stream, name, histogram);
        }

        return;
    }

    /**
     * Print MultivariateTracker matrix values and their dimension titles
     * @param[out] stream Destination stream
     * @param[in] tracker Source tracker
     * @param[in] title Matrix title
     * @param[in] matrix Pointer to matrix function
     */
    template<class T>
    void
    writeMatrix(std::ostream & stream,
                MultivariateTracker & tracker,
                std::string const & title,
                T (MultivariateTracker::*matrix)(long const i, long const j) const)
    {

        int numDim = tracker.getDimension();

        // Compute column widths
        std::vector<long> width(numDim,0);
        for (long i = 0; i < numDim; i++)
        {
            StringParser::updateWidthString(tracker.getName(i), width[i]);
        }
        for (long i = 0; i < numDim; i++)
        {
            for (long j = 0; j < numDim; j++)
            {
                StringParser::updateWidthNumber<T>((tracker.*matrix)(i,j), width[j]);
            }
        }

        long totalWidth = 0;
        for (std::vector<long>::iterator it = width.begin(); it != width.end(); ++it)
        {
            (*it) += 3;
            totalWidth += *it;
        }
        totalWidth = std::max<long>(totalWidth, 20);

        // Print data headers
        std::string header(totalWidth,'='); 
        std::string footer(totalWidth,'-'); 
        stream << this->_comment << header << std::endl;
        stream << this->_comment << std::setfill (' ') << std::setw (totalWidth/2) << title << std::endl;
        stream << this->_comment << header << std::endl;
        stream << this->_comment;
        if (numDim == 0)
        {
            stream << "No Data" << std::endl;
            return;
        }
        for (long i = 0; i < numDim; i++)
        {
            stream << std::right << std::setw(width[i]) << tracker.getName(i);
        }
        stream << std::endl;
        stream << this->_comment << footer << std::endl;

        // Print data
        std::string indent(this->_comment.size(),' '); 
        for (long i = 0; i < numDim; i++)
        {
            stream << indent;
            for (long j = 0; j < numDim; j++)
            {
                stream << std::right << std::setw(width[j]) << StringParser::parseNumber<T>((tracker.*matrix)(i,j));
            }
            stream << std::endl;
        }
        stream << std::endl;
    
        return;
        
    }

    /**
     * Print summary statistics
     * @param[out] stream Destination stream
     * @param[in] tracker Source tracker
     */
    void
    writeStatistics(std::ostream & stream,
                    MultivariateTracker & tracker)
    {

        int numDim = tracker.getDimension();
    
        // Compute column widths
        std::vector<long> width(6,0);
        StringParser::updateWidthString("Dimension", width[0]);
        StringParser::updateWidthString("Count", width[1]);
        StringParser::updateWidthString("Minimum", width[2]);
        StringParser::updateWidthString("Mean", width[3]);
        StringParser::updateWidthString("Maximum", width[4]);
        StringParser::updateWidthString("Stdev", width[5]);
        for (long i = 0; i < numDim; i++)
        {
            StringParser::updateWidthString(tracker.getName(i), width[0]);
            StringParser::updateWidthNumber<long>(tracker.getCount(i), width[1]);
            StringParser::updateWidthString(StringParser::parseStatistic<double>(tracker.getCount(i), tracker.getMinimum(i)), width[2]);
            StringParser::updateWidthString(StringParser::parseStatistic<double>(tracker.getCount(i), tracker.getMean(i)), width[3]);
            StringParser::updateWidthString(StringParser::parseStatistic<double>(tracker.getCount(i), tracker.getMaximum(i)), width[4]);
            StringParser::updateWidthString(StringParser::parseStatistic<double>(tracker.getCount(i), sqrt(tracker.getVariance(i))), width[5]);
        }

        long totalWidth = 0;
        for (std::vector<long>::iterator it = width.begin(); it != width.end(); ++it)
        {
            (*it) += 3;
            totalWidth += *it;
        }

        // Print data headers
        std::string header(totalWidth,'='); 
        std::string footer(totalWidth,'-'); 
        stream << this->_comment << header << std::endl;
        stream << this->_comment << std::setfill (' ') << std::setw (totalWidth/2) << "Statistics" << std::endl;
        stream << this->_comment << header << std::endl;
        stream << this->_comment;
        stream << std::right << std::setw(width[0]) << "Dimension";
        stream << std::right << std::setw(width[1]) << "Count";
        stream << std::right << std::setw(width[2]) << "Minimum";
        stream << std::right << std::setw(width[3]) << "Mean";
        stream << std::right << std::setw(width[4]) << "Maximum";
        stream << std::right << std::setw(width[5]) << "Stdev";
        stream << std::endl;
        stream << this->_comment << footer << std::endl;
        
        // Print data
        std::string indent(this->_comment.size(),' '); 
        for (long i = 0; i < numDim; i++)
        {
            stream << indent;
            stream << std::right << std::setw(width[0]) << tracker.getName(i);
            stream << std::right << std::setw(width[1]) << StringParser::parseNumber<long>(tracker.getCount(i));
            stream << std::right << std::setw(width[2]) << StringParser::parseStatistic<double>(tracker.getCount(i), tracker.getMinimum(i));
            stream << std::right << std::setw(width[3]) << StringParser::parseStatistic<double>(tracker.getCount(i), tracker.getMean(i));
            stream << std::right << std::setw(width[4]) << StringParser::parseStatistic<double>(tracker.getCount(i), tracker.getMaximum(i));
            stream << std::right << std::setw(width[5]) << StringParser::parseStatistic<double>(tracker.getCount(i), sqrt(tracker.getVariance(i)));
            stream << std::endl;
        }
        stream << std::endl;

        return;
        
    }
    
    /**
     * Write histogram to a stream
     * @param[out] stream Destination stream
     * @param[in] tracker Source tracker
     * @param[in] histogram Histogram to write
     */
    template <class T>
    void
    writeHistogram(std::ostream & stream,
                   std::string const & title,
                   DynamicHistogram<T> & histogram)
    {

        try
        {

            unsigned int nBins = 0;
            std::vector<double> freq;
            std::vector<double> pdf;
            std::vector<double> cdf;

            // Compute column widths
            std::vector<long> width(6,0);
            StringParser::updateWidthString("Bin Min", width[0]);
            StringParser::updateWidthString("Bin Max", width[1]);
            StringParser::updateWidthString("Center", width[2]);
            StringParser::updateWidthString("Frequency", width[3]);
            StringParser::updateWidthString("PDF", width[4]);
            StringParser::updateWidthString("CDF", width[5]);

            nBins = histogram.count();
            freq = histogram.frequencies();
            pdf = histogram.pdf();
            cdf = histogram.cdf();
            for (unsigned int i = 0; i < nBins; i++)
            {
                StringParser::updateWidthNumber<long>(i, width[0]);
                StringParser::updateWidthNumber<double>(histogram.bin(i)-histogram.width(i)/2, width[1]);
                StringParser::updateWidthNumber<double>(histogram.bin(i)+histogram.width(i)/2, width[2]);
                StringParser::updateWidthNumber<long long>((long long)freq.at(i), width[3]);
                StringParser::updateWidthNumber<double>(pdf.at(i), width[4]);
                StringParser::updateWidthNumber<double>(cdf.at(i), width[5]);
            }

            long totalWidth = 0;
            for (std::vector<long>::iterator it = width.begin(); it != width.end(); ++it)
            {
                (*it) += 3;
                totalWidth += *it;
            }

            // Print data headers
            std::string header(totalWidth,'=');
            std::string footer(totalWidth,'-');
            stream << this->_comment << header << std::endl;
            stream << this->_comment << std::setfill (' ') << std::setw (totalWidth/2) << "Histogram of " << title << std::endl;
            stream << this->_comment << header << std::endl;
            stream << this->_comment;
            stream << std::right << std::setw(width[0]) << "Bin"
                   << std::right << std::setw(width[1]) << "Min"
                   << std::right << std::setw(width[2]) << "Max"
                   << std::right << std::setw(width[3]) << "Frequency"
                   << std::right << std::setw(width[4]) << "PDF"
                   << std::right << std::setw(width[5]) << "CDF"
                   << std::endl;

            // Print data
            std::string indent(this->_comment.size(), ' ');
            for (unsigned int i = 0; i < nBins; i++)
            {
                stream << indent;
                stream << std::right << std::setw(width[0]) << StringParser::parseNumber<long>(i);
                stream << std::right << std::setw(width[1]) << StringParser::parseNumber<double>(histogram.bin(i)-histogram.width(i)/2);
                stream << std::right << std::setw(width[2]) << StringParser::parseNumber<double>(histogram.bin(i)+histogram.width(i)/2);
                stream << std::right << std::setw(width[3]) << StringParser::parseNumber<long long>((long long)freq.at(i));
                stream << std::right << std::setw(width[4]) << StringParser::parseNumber<double>(pdf.at(i));
                stream << std::right << std::setw(width[5]) << StringParser::parseNumber<double>(cdf.at(i));
                stream << std::endl;
            }

        }
        catch (...)
        {
            stream << "No valid data" << std::endl;
        }
        stream << std::endl;

        return;
    }

private:

    std::string _comment;
    
};

/**
 * @class RowSampler
 * @brief Class to decide if current entry should be sampled
 */
class RowSampler
{
public:

    /**
     * @param[in] options Sampler Options
     */
    RowSampler(SamplerOptions const & options) :
        _count(0),
        _next(1),
        _sampled(0),
        _mode(options.mode),
        _step(options.step)
    {
    }

    /**
     * Determine if current sample should be used
     * @return True if sample should be used otherwise false
     */
    bool
    sample()
    {
        // Don't bother if samling every row
        if (this->_step == 1) return true;

        // Bump row counter
        this->_count++;

        // Check if this row is sampled
        if (this->_count == this->_next)
        {
            // Bump sampled counter
            this->_sampled++;

            // If yes, determine next row in this sampling bin to sample
            unsigned int row = 1;
            if (this->_mode == SamplerOptions::RANDOM)
            {
                row = 1 + (unsigned int)(rand() % this->_step);
            }
            this->_next = this->_sampled * this->_step + row;

            return true;
        }

        return false;

    }

private:

     /**
      * Number of rows tested for sampling
      */
    unsigned int _count;

     /**
      * Next row to sample
      */
    unsigned int _next;

   /**
    * Sampling mode
    */
    SamplerOptions::SampleMode _mode;

    /**
     * Number of sampled rows
     */
    unsigned int _sampled;

    /**
     * Sampling step size
     */
    unsigned int _step;

};

/**
 * @class StatisticsApp
 * @brief statistics application
 */
class StatisticsApp
{

public:

    /**
     * Create application
     */
    StatisticsApp(CommandLineOptions const & options) :
        _options(options),
        _outputStream(0),
        _tracker(MultivariateTracker(0, options.statisticsOptions))
    {
        Logger::logLevel = (Logger::Level::Type)this->_options.verboseLevel;

        // ==================================================
        // Seed random number generator
        // --------------------------------------------------
        srand((options.seed == 0 ? (unsigned int)time(0) : options.seed));

        // ==================================================
        // Define output stream
        // --------------------------------------------------
        this->_outputFile.open(this->_options.fileOutput.c_str(), std::ios::out);
        if (!this->_outputFile.good() && !this->_options.fileOutput.empty())
        {
            LOG_MESSAGE(Logger::Level::ERROR, "Output file failed to open. Writting to standard output.");
        }
        this->_outputStream = (this->_outputFile.good() ? &this->_outputFile : &std::cout);
    }

    /**
     * Run application
     * @return Status code
     */
    int
    run(void)
    {

        // ==================================================
        // Application parameters
        // --------------------------------------------------
        int numDim = 0;
        int numDimMask = 0;

        // ==================================================
        // Define input stream
        // --------------------------------------------------
        std::ifstream fileInput(this->_options.fileInput.c_str(), std::ios::in);
        std::istream * source = (fileInput.good() ? &fileInput : &std::cin);

        LOG_MESSAGE(Logger::Level::INFO, ApplicationProperties::NAME << " is starting");

        // ==================================================
        // Define output writer
        // --------------------------------------------------
        StatisticsWriter writer;

        // ==================================================
        // Define row sampler
        // --------------------------------------------------
        RowSampler sampler(this->_options.sampling);

        // ==================================================
        // Loop through source data
        // --------------------------------------------------
        
        // titles is populated by parsed line if user provides a titles row, otherwise indices are used
        std::vector<std::string> titles;
        StringSplitter parser;
        
        // Counters used to keep track of line reads if reshaping is requested by user
        unsigned int counter = 0;
        unsigned int counterSinceLastParsing = 0;
        std::string lineToReshape = "";
        
        // DataVector contains the parsed data from each line.
        // It's size will get set on the first successful parsing.
        DataVector data;
        DataVector dataMasked;
        
        bool reading = true;
        while (reading)
        {

            // ==================================================
            // Extract current line
            // --------------------------------------------------
            std::string line = StatisticsApp::readline(*source);
            if (!(*source).good())
            {
                LOG_MESSAGE(Logger::Level::DETAIL, "End of file found. Exiting.");
                break;
            }

            // ==================================================
            // Check if line is empty and if we need to exit
            // --------------------------------------------------
            if (line.size() == 0)
            {
                if (this->_options.blankEOF)
                {
                    LOG_MESSAGE(Logger::Level::DETAIL, "Blank line found. Exiting.");
                    break;
                }
                else
                {
                    LOG_MESSAGE(Logger::Level::DETAIL, "Blank line found. Skipping.");
                    continue;
                }
            }
            
            // ==================================================
            // Check if row is comment and needs to be skipped
            // --------------------------------------------------
            if (line.compare(0,1,this->_options.comment.c_str()) == 0) continue;
            
            // ==================================================
            // Check if row needs to be reshaped into column
            // --------------------------------------------------

            // If yes, then append current line to previous
            // lines with delimiter and continue to next read
            counterSinceLastParsing++;
            lineToReshape += (counterSinceLastParsing == 1 ? "" : this->_options.delimiter) + line;
            if (counterSinceLastParsing < this->_options.numLinesToReshape) continue;

            // Prepare line and counter for parsing
            line = lineToReshape;
            counter++;

            // Reset counter and appended line for next read
            counterSinceLastParsing = 0;
            lineToReshape = "";

            // ==================================================
            // Extract header titles
            // --------------------------------------------------
            if (this->_options.headerRow == counter)
            {
                parser.initialize(line, this->_options.delimiter, this->_options.removeDuplicates, this->_options.filterColumns);
                titles = parser.tokens();
                continue;
            }
            
            // ==================================================
            // Skip header rows
            // --------------------------------------------------
            if (counter <= this->_options.numLinesToSkip)
            {
                LOG_MESSAGE(Logger::Level::INFO, "Skipping header line.");
                continue;
            }

            // ==================================================
            // Skip unneeded rows
            // --------------------------------------------------
            if (counter > (this->_options.numLinesToSkip + this->_options.numLinesToKeep))
            {
                LOG_MESSAGE(Logger::Level::INFO, "Skipping remaining lines.");
                break;
            }

            // ==================================================
            // Skip unsampled rows
            // --------------------------------------------------
            if (!sampler.sample())
            {
                LOG_MESSAGE(Logger::Level::DETAIL, "Row not sampled. Skipping");
                continue;
            }

            // ==================================================
            // Parse current line
            // --------------------------------------------------
            parser.initialize(line, this->_options.delimiter, this->_options.removeDuplicates);
            StringSplitter::size_type numTokens = parser.size();
            LOG_MESSAGE(Logger::Level::DETAIL, "Line = " << counter << " - content =\"" << line << "\"");

            // ==================================================
            // Check and compute source data dimensions
            // --------------------------------------------------
            bool isDimSet = numDim > 0;
            if (!isDimSet)
            {
                numDim = (int) numTokens;
                data.resize(numDim);

                // Make sure we only keep column filters if they are within actual number of columns
                this->_options.filterColumns.resize(numDim, (this->_options.filterColumns.size() == 0));

                // Compute how many columns are masked out
                numDimMask = std::accumulate(this->_options.filterColumns.begin(), this->_options.filterColumns.end(), 0);
                dataMasked.resize(numDimMask);

                LOG_MESSAGE(Logger::Level::INFO, "Data dimensions set to " << numDim);
                LOG_MESSAGE(Logger::Level::INFO, "Mask dimensions set to " << numDimMask);
            }

            bool isDimFailed = (numDim != (int) numTokens);
            if (isDimFailed)
            {
                LOG_MESSAGE(Logger::Level::WARNING, "Line = " << counter << " - Statistics dimensions does not equal data dimensions." << (this->_options.strictParsing ? " Skipping." : ""));
                if (this->_options.strictParsing)
                {
                    continue;
                }
            }

            LOG_MESSAGE(Logger::Level::DEBUG, "Line = " << counter << " - Number of tokens found: " << numTokens);

            // ==================================================
            // Apply string filters
            // --------------------------------------------------
            bool isStringFiltered = this->_options.filterRows.isFiltered(parser.get());
            if (!isStringFiltered)
            {
                LOG_MESSAGE(Logger::Level::DEBUG, "Line = " << counter << " - filtered out by string filter");
                continue;
            }

            // ==================================================
            // Convert tokens into data array
            // --------------------------------------------------
            for (StringSplitter::size_type i = 0; i < (StringSplitter::size_type) numDim; ++i)
            {

                if (numTokens <= i) break;

                data.at(i).active = parser.toValue<double>(i, data.at(i).value);

                LOG_MESSAGE(Logger::Level::DETAIL, "Line = " << counter << " - " << (data.at(i).active ? "Valid token" : "Invalid token") << " = " << parser.at(i));

            }

            // ==================================================
            // Apply numeric filters
            // --------------------------------------------------
            bool isNumericFiltered = this->_options.filterRows.isFiltered(data);
            if (!isNumericFiltered)
            {
                LOG_MESSAGE(Logger::Level::DEBUG, "Line = " << counter << " - filtered out by numeric filter");
                continue;
            }

            // ==================================================
            // Apply column filters
            // --------------------------------------------------
            DataVector::copy(data, numDim, dataMasked, numDimMask, this->_options.filterColumns, numDim);

            // ==================================================
            // Write current data
            // --------------------------------------------------
            if (this->_options.showFilteredData)
            {
                writer.writeData(*this->_outputStream, dataMasked);
            }

            // ==================================================
            // Update statistics with current data
            // --------------------------------------------------
            this->_tracker.update(dataMasked);

            // Since the data vector is preallocated, we deactivate all entries in case they don't exist in the next parsed line
            // This way they don't get inadvertently included in the next line's statistics
            data.deactivate();

        }

        // ==================================================
        // Update statistics names
        // --------------------------------------------------
        if (numDimMask != (long)titles.size())
        {
            // Show error message if user provided a bad header row
            if (this->_options.headerRow > 0)
            {
                LOG_MESSAGE(Logger::Level::ERROR, "Size of header titles does not match data dimensions so using indices as header titles.");
            }
            titles.resize(numDimMask,"");
            for (long i = 0; i < numDimMask; i++)
            {
                titles.at(i) = StringParser::parseNumber<long>(i+1);
            }
        }
        for (long i = 0; i < numDimMask && i < this->_tracker.getDimension(); i++)
        {
            this->_tracker.setName(i, titles.at(i));
        }

        LOG_MESSAGE(Logger::Level::INFO, ApplicationProperties::NAME << " is complete");
    
        return 0;
        
    }

    /**
     * Dispaly output statistics
     * @return True if statistics displayed successfully otherwise false
     */
    bool
    display()
    {

        // ==================================================
        // Write output to stream
        // --------------------------------------------------
        StatisticsWriter writer;
        
        // Show statistics
        if (this->_options.showStatistics)
        {
            writer.writeStatistics(*this->_outputStream, this->_tracker);
        }
        
        // Show covariance
        if (this->_options.showCovariance)
        {
            writer.writeMatrix<double>(*this->_outputStream, this->_tracker, "Covariance", &MultivariateTracker::getCovariance);
        }

        // Show correlation
        if (this->_options.showCorrelation)
        {
            writer.writeMatrix<double>(*this->_outputStream, this->_tracker, "Correlation", &MultivariateTracker::getCorrelation);
        }

        // Show least squares fit offset
        if (this->_options.showLeastSquaresOffset)
        {
            writer.writeMatrix<double>(*this->_outputStream, this->_tracker, "Offset", &MultivariateTracker::getLeastSquaresOffset);
        }

        // Show least squares fit slope
        if (this->_options.showLeastSquaresSlope)
        {
            writer.writeMatrix<double>(*this->_outputStream, this->_tracker, "Slope", &MultivariateTracker::getLeastSquaresSlope);
        }
        
        // Show histogram
        if (this->_options.showHistogram)
        {
            writer.writeHistograms(*this->_outputStream, this->_tracker);
        }

        return true;

    }

private:

    /**
     * @brief Read current line from stream
     * @param[in] stream Source stream
     * @return Read line as string
     */
    static
    std::string
    readline(std::istream & stream)
    {
    
        // Output line
        std::string line;

        // Get stream buffer and guard with sentry
        std::istream::sentry guard(stream, true);
        std::streambuf * buffer = stream.rdbuf();
        
        // Iterate over source buffer, append new characters to
        // destination string, and check for the following
        // line endings:
        //  - LF = Linux
        //  - CR = Apple
        //  - CRLF = Windows
        //  - EOF = End-Of-file
        while (true)
        {
        
            // Extract current character and bump stream pointer
            int character = buffer->sbumpc();
            
            // Check if current character is a line feed
            if (character ==  StatisticsApp::LF)
            {
                return line;
            }
            
            // Check if current character is a carriage return
            if (character == StatisticsApp::CR)
            {
                // Get current character without bumping stream pointer
                character = buffer->sgetc();
                if(character == StatisticsApp::LF)
                {
                    // Remove line feed from stream and bump stream pointer
                    buffer->sbumpc();
                }
                return line;
            }
            
            // Check if current character is the end of file
            if (character == StatisticsApp::EF)
            {
                if(line.empty())
                {
                    stream.setstate(std::ios::eofbit);
                }
                return line;
            }
            
            // Append current character to destination string
            line += (char) character;
            
        }
        
        return line;
        
    }
    
private:

    /**
     * Line feed character
     */
    static const int LF = '\n';
    
    /**
     * Carriage return character
     */
    static const int CR = '\r';

    /**
     * End Of File character
     */
    static const int EF = EOF;

private:

    /**
     * Options data delimiter
     */
    CommandLineOptions _options;

    /**
     * Output file
     */
    std::ofstream _outputFile;

    /**
     * Output stream
     */
    std::ostream * _outputStream;

    /**
     * Statistics
     */
    MultivariateTracker _tracker;

};

/**
 * Application status codes
 */
namespace AppStatus
{
    enum Status
    {
        SUCCESS = 0,
        FAILED_PARSING,
        FAILED_RUN,
        FAILED_DISPLAY
    };
}

/**
 * @brief program driver
 * @param[in] argc Number of arguments including executable path
 * @param[in] argv Character array of arguments
 * @return Exit status
 *          - 0 = Success
 *          - 1 = Argument parsing failed
 *          - 2 = Computing statistics failed
 *          - 3 = Displaying statistics failed
 */
int
pdrv(int argc,
     char * argv[])
{

    // ======================================================================
    // Parse command line arguments
    // ----------------------------------------------------------------------
    CommandLineOptions options;

    try
    {
        CommandLineParser parser(argc, argv);
        if (parser.showUsage())
        {
            CommandLineParser::printUsage();
            return AppStatus::SUCCESS;
        }
        else if (parser.showVersion())
        {
            CommandLineParser::printVersion();
            return AppStatus::SUCCESS;
        }
        options = parser.options;
    }
    catch (std::exception & ex)
    {
        LOG_MESSAGE(Logger::Level::FATAL, ex.what());
        return AppStatus::FAILED_PARSING;
    }

    // ======================================================================
    // Launch application
    // ----------------------------------------------------------------------
    StatisticsApp application(options);

    try
    {
        int status = application.run();
        if (status != 0) return AppStatus::FAILED_RUN;
    }
    catch (std::exception & ex)
    {
        LOG_MESSAGE(Logger::Level::FATAL, ex.what());
        return AppStatus::FAILED_RUN;
    }

    // ======================================================================
    // Display results
    // ----------------------------------------------------------------------
    
    try
    {
        bool isDisplayValid = application.display();
        if (!isDisplayValid) return AppStatus::FAILED_DISPLAY;
    }
    catch (std::exception & ex)
    {
        LOG_MESSAGE(Logger::Level::FATAL, ex.what());
        return AppStatus::FAILED_DISPLAY;
    }
        
    return AppStatus::SUCCESS;

}

#ifndef _CLISTATS_TESTING
int
main(int argc,
     char * argv[])
{
    return pdrv(argc, argv);
}
#endif
