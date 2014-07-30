
/**
 * Unit tests for clistats
 * @author Daniel Pulido <dpmcmlxxvi@gmail.com>
 * @copyright Copyright (c) 2014 Daniel Pulido <dpmcmlxxvi@gmail.com>
 * @file tests.cpp
 * @license MIT License (http://opensource.org/licenses/MIT)
 */
 
#define _CLISTATS_TESTING

#include "../src/clistats.cpp"

// macros for asserting boolean condition with messae
#define ASSERT(condition, message)                  \
    {                                               \
        if (!(condition))                           \
        {                                           \
            std::stringstream msg;                  \
            msg << message << std::endl;            \
            throw std::runtime_error(msg.str());    \
        }                                           \
    }
    
#define ASSERT_EQUAL(value1,value2, message)        \
    {                                               \
        ASSERT((value1)==(value2), message)         \
    }

#define ASSERT_EXCEPTION(condition, message)        \
    {                                               \
        try                                         \
        {                                           \
            condition;                              \
            ASSERT(false, message);                 \
        }                                           \
        catch (...)                                 \
        {                                           \
        }                                           \
    }

/**
 * Test StringParser::parseNumber
 */
void test_stringparser_parsenumer()
{

    ASSERT_EQUAL(StringParser::parseNumber<short>(1), "1", "StringParser::parseNumber<short> failed on short");
    ASSERT_EQUAL(StringParser::parseNumber<int>(1), "1", "StringParser::parseNumber<int> failed on parsing int");
    ASSERT_EQUAL(StringParser::parseNumber<long>(1), "1", "StringParser::parseNumber<long> failed on parsing long");
    ASSERT_EQUAL(StringParser::parseNumber<long long>(1), "1", "StringParser::parseNumber<long long> failed on parsing long");
    ASSERT_EQUAL(StringParser::parseNumber<float>((float)1.1), "1.100000", "StringParser::parseNumber<float> failed on parsing float");
    ASSERT_EQUAL(StringParser::parseNumber<double>(1.1), "1.100000", "StringParser::parseNumber<double> failed on parsing double");
    ASSERT_EQUAL(StringParser::parseNumber<long double>(1.1), "1.100000", "StringParser::parseNumber<long double> failed on parsing long double");
    ASSERT_EQUAL(StringParser::parseNumber<float>(std::numeric_limits<float>::quiet_NaN()), "nan", "StringParser::parseNumber<float> failed on parsing nan");
    ASSERT_EQUAL(StringParser::parseNumber<double>(std::numeric_limits<double>::quiet_NaN()), "nan", "StringParser::parseNumber<double> failed on parsing nan");

}

/**
 * Test StringParser::parseStatistic
 */
void test_stringparser_parsestatistic()
{

    ASSERT_EQUAL(StringParser::parseStatistic<short>(0,1), "nan", "StringParser::parseStatistic<short> failed on short");
    ASSERT_EQUAL(StringParser::parseStatistic<int>(0,1), "nan", "StringParser::parseStatistic<int> failed on parsing int");
    ASSERT_EQUAL(StringParser::parseStatistic<long>(0,1), "nan", "StringParser::parseStatistic<long> failed on parsing long");
    ASSERT_EQUAL(StringParser::parseStatistic<long long>(0,1), "nan", "StringParser::parseStatistic<long long> failed on parsing long");
    ASSERT_EQUAL(StringParser::parseStatistic<float>(0,(float)1.1), "nan", "StringParser::parseStatistic<float> failed on parsing float");
    ASSERT_EQUAL(StringParser::parseStatistic<double>(0,1.1), "nan", "StringParser::parseStatistic<double> failed on parsing double");
    ASSERT_EQUAL(StringParser::parseStatistic<long double>(0,1.1), "nan", "StringParser::parseStatistic<long double> failed on parsing long double");

}

/**
 * Test StringParser::replacen
 */
void test_stringparser_replacen()
{
    std::string content = "barbar";
    StringParser::replacen(content,"bar","foo");
    ASSERT_EQUAL(content, "foobar", "StringParser::replacen failed barbar -> foobar");
    StringParser::replacen(content,"bar","foo",1);
    ASSERT_EQUAL(content, "foofoo", "StringParser::replacen failed foobar -> foofoo");
    StringParser::replacen(content,"foo","bar",2);
    ASSERT_EQUAL(content, "barbar", "StringParser::replacen failed foofoo -> barbar");

}

/**
 * Test StringParser::toValue
 */
void test_stringparser_tovalue()
{

    short shortValue;
    ASSERT(!StringParser::toValue<short>("1.1", shortValue), "StringParser::toValue failed on parsing invalid short");
    ASSERT(StringParser::toValue<short>("1", shortValue), "StringParser::toValue failed on parsing short");
    ASSERT_EQUAL(shortValue, 1, "StringParser::toValue failed on converting string to short value");

    int intValue;
    ASSERT(!StringParser::toValue<int>("1.1", intValue), "StringParser::toValue failed on parsing invalid int");
    ASSERT(StringParser::toValue<int>("1", intValue), "StringParser::toValue failed on parsing int");
    ASSERT_EQUAL(intValue, 1, "StringParser::toValue failed on converting string to int value");

    long longValue;
    ASSERT(!StringParser::toValue<long>("1.1", longValue), "StringParser::toValue failed on parsing invalid long");
    ASSERT(StringParser::toValue<long>("1", longValue), "StringParser::toValue failed on parsing long");
    ASSERT_EQUAL(longValue, 1, "StringParser::toValue failed on converting string to long value");

    long long longlongValue;
    ASSERT(!StringParser::toValue<long long>("1.1", longlongValue), "StringParser::toValue failed on parsing invalid long long");
    ASSERT(StringParser::toValue<long long>("1", longlongValue), "StringParser::toValue failed on parsing long long");
    ASSERT_EQUAL(longlongValue, 1, "StringParser::toValue failed on converting string to long long value");

    float floatValue;
    ASSERT(!StringParser::toValue<float>("a", floatValue), "StringParser::toValue failed on parsing invalid float");
    ASSERT(StringParser::toValue<float>("1.1", floatValue), "StringParser::toValue failed on parsing float");
    ASSERT_EQUAL(floatValue, (float)1.1, "StringParser::toValue failed on converting string to float value");

    double doubleValue;
    ASSERT(!StringParser::toValue<double>("a", doubleValue), "StringParser::toValue failed on parsing invalid double");
    ASSERT(StringParser::toValue<double>("1.1", doubleValue), "StringParser::toValue failed on parsing double");
    ASSERT_EQUAL(doubleValue, (double)1.1, "StringParser::toValue failed on converting string to double value");

}

/**
 * Test StringParser::trim
 */
void test_stringparser_trim()
{

    ASSERT_EQUAL(StringParser::trim("    foo bar"), "foo bar", "StringParser::trim failed \"    foo bar\" -> \"foo bar\"");
    ASSERT_EQUAL(StringParser::trim("foo bar    "), "foo bar", "StringParser::trim failed \"foo bar    \" -> \"foo bar\"");
    ASSERT_EQUAL(StringParser::trim("    foo bar    "), "foo bar", "StringParser::trim failed \"    foo bar    \" -> \"foo bar\"");

}

/**
 * Test StringParser::toValue
 */
void test_stringparser_updatewidthnumber()
{

    long width = 0;
    StringParser::updateWidthNumber<short>(1, width);
    ASSERT_EQUAL(width, 1, "StringParser::updateWidthNumber failed on updating short");
    StringParser::updateWidthNumber<int>(1, width);
    ASSERT_EQUAL(width, 1, "StringParser::updateWidthNumber failed on updating int");
    StringParser::updateWidthNumber<long>(1, width);
    ASSERT_EQUAL(width, 1, "StringParser::updateWidthNumber failed on updating long");
    StringParser::updateWidthNumber<long long>(1, width);
    ASSERT_EQUAL(width, 1, "StringParser::updateWidthNumber failed on updating long long");
    StringParser::updateWidthNumber<float>(1, width);
    ASSERT_EQUAL(width, 8, "StringParser::updateWidthNumber failed on updating float");
    StringParser::updateWidthNumber<double>(1, width);
    ASSERT_EQUAL(width, 8, "StringParser::updateWidthNumber failed on updating double");

}

/**
 * Test StringParser::toValue
 */
void test_stringparser_updatewidthstring()
{

    long width = 0;
    StringParser::updateWidthString("1", width);
    ASSERT_EQUAL(width, 1, "StringParser::updateWidthNumber failed on updating string \"1\"");
    StringParser::updateWidthString("12", width);
    ASSERT_EQUAL(width, 2, "StringParser::updateWidthNumber failed on updating string \"12\"");
    StringParser::updateWidthString("", width);
    ASSERT_EQUAL(width, 2, "StringParser::updateWidthNumber failed on updating string \"\"");
    StringParser::updateWidthString("123", width);
    ASSERT_EQUAL(width, 3, "StringParser::updateWidthNumber failed on updating string \"123\"");

}

/**
 * Test StringSplitter
 */
void test_stringsplitter()
{

    int count = 3;
    std::string valuesstring = "foo,bar,1";
    const char * valueschar[] = {"foo", "bar", "1"};

    StringSplitter splitter(valuesstring,",");
    for (int i = 0; i < count; i++)
    {
        ASSERT_EQUAL(splitter.at(i), valueschar[i], "StringSplitter::at failed on " << valueschar[i]);
    }

    std::vector<std::string> valuesGet(splitter.get());
    for (int i = 0; i < count; i++)
    {
        ASSERT_EQUAL(valueschar[i], valuesGet.at(i), "StringSplitter::get failed at " << valueschar[i]);
    }

    ASSERT_EQUAL(splitter.size(), count, "StringSplitter::size failed");

    std::vector<std::string> valuesTokens(splitter.tokens());
    for (int i = 0; i < count; i++)
    {
        ASSERT_EQUAL(valueschar[i], valuesTokens.at(i), "StringSplitter::tokens failed at " << valueschar[i]);
    }

    int value;
    ASSERT_EXCEPTION(splitter.toValue(0,value), "StringSplitter::toValue failed to not parse non-numeric value");
    ASSERT_EXCEPTION(splitter.toValue(1,value), "StringSplitter::toValue failed to not parse non-numeric value");
    ASSERT(splitter.toValue(2,value), "StringSplitter::toValue failed to parse numeric value");

}

/**
 * Test StringSplitter::toIntegers
 */
void test_stringsplitter_tointegers()
{

    std::vector<int> integersTest;
    ASSERT(!StringSplitter::toIntegers("1:",integersTest), "StringSplitter::toIntegers failed to not parse invalid range string");
    ASSERT(StringSplitter::toIntegers("1:3",integersTest), "StringSplitter::toIntegers failed to parse range string");
    ASSERT_EQUAL(integersTest.size(), 3, "StringSplitter::toIntegers failed to convert range string to correct vector size");
    ASSERT_EQUAL(integersTest.at(0), 1, "StringSplitter::toIntegers failed to convert range string value 1");
    ASSERT_EQUAL(integersTest.at(1), 2, "StringSplitter::toIntegers failed to convert range string value 2");
    ASSERT_EQUAL(integersTest.at(2), 3, "StringSplitter::toIntegers failed to convert range string value 3");
    
}

/**
 * Test FixedSizeCache
 */
void test_fixedsizecache()
{

    int size = 3;
    FixedSizeCache<int> cache(size);
    
    ASSERT(cache.empty(), "FixedSizeCache::empty failed on empty cache");
    ASSERT(!cache.full(), "FixedSizeCache::full failed on non-full cache");
    ASSERT_EXCEPTION(cache.max(), "FixedSizeCache::max failed on empty cache");
    ASSERT_EXCEPTION(cache.min(), "FixedSizeCache::min failed on empty cache");

    cache.add(1);
    cache.add(2);
    cache.add(3);
    ASSERT_EXCEPTION(cache.add(4), "FixedSizeCache::add failed to not add on full cache");
    for (int i = 0; i < size; i++)
    {
        ASSERT_EQUAL(cache.at(i), i+1, "FixedSizeCache::at failed on " << i);
    }
    ASSERT_EXCEPTION(cache.at(3), "FixedSizeCache::at failed to not access on full cache");
    ASSERT_EQUAL(cache.count(), size, "FixedSizeCache::count failed");
    ASSERT(!cache.empty(), "FixedSizeCache::count failed on non-empty cache");
    ASSERT(cache.full(), "FixedSizeCache::full failed on full cache");
    
    cache.reset();
    ASSERT(cache.empty(), "FixedSizeCache::empty failed on reset cache");
    ASSERT(!cache.full(), "FixedSizeCache::full failed on reset cache");
    ASSERT_EXCEPTION(cache.max(), "FixedSizeCache::max failed on reset cache");
    ASSERT_EXCEPTION(cache.min(), "FixedSizeCache::min failed on reset cache");
    
}

/**
 * Test DynamicHistogram disabled
 */
void test_dynamichistogram_disabled()
{

    DynamicHistogramOptions options(false);
    DynamicHistogram<double> histogramDisabled(options);
    ASSERT(!histogramDisabled.add(1), "DynamicHistogram::add failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.bin(0), "DynamicHistogram::bin failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.cdf(), "DynamicHistogram::cdf failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.contains(0), "DynamicHistogram::contains failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.count(), "DynamicHistogram::count failed on disabled histogram");
    ASSERT(!histogramDisabled.enabled(), "DynamicHistogram::enabled failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.frequency(0), "DynamicHistogram::frequency failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.frequencies(), "DynamicHistogram::frequencies failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.index(0), "DynamicHistogram::index failed on disabled histogram");
    ASSERT(!histogramDisabled.initialized(), "DynamicHistogram::initialized failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.max(), "DynamicHistogram::max failed on disabled histogram");
    ASSERT_EQUAL(histogramDisabled.merges(), 0, "DynamicHistogram::merges failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.min(), "DynamicHistogram::min failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.order(0), "DynamicHistogram::order failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.pdf(), "DynamicHistogram::pdf failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.total(), "DynamicHistogram::total failed on disabled histogram");
    ASSERT_EXCEPTION(histogramDisabled.width(0), "DynamicHistogram::width failed on disabled histogram");

}

/**
 * Test DynamicHistogram enabled
 */
void test_dynamichistogram_enabled()
{

    DynamicHistogramOptions options(false);
    DynamicHistogram<double> histogramDisabled(options);
    options.enabled = true;
    options.binCount = 4;
    options.cacheSize = 4;
    const double valuesdouble[] = {0.0,1.0,2.0,4.0};
    
    DynamicHistogram<double> histogramEnabled(options);

    for (int i = 0; i < options.binCount; i++)
    {
        ASSERT(!histogramEnabled.add(valuesdouble[i]), "DynamicHistogram::add failed on adding to cache value = " << i);
    }
    for (int i = 0; i < options.binCount; i++)
    {
        ASSERT(histogramEnabled.add(valuesdouble[i]), "DynamicHistogram::add failed on adding to histogram value = " << i);
    }

    ASSERT_EQUAL(histogramEnabled.bin(0), 0.5, "DynamicHistogram::bin failed on enabled histogram at 0");
    ASSERT_EQUAL(histogramEnabled.bin(1), 1.5, "DynamicHistogram::bin failed on enabled histogram at 1");
    ASSERT_EQUAL(histogramEnabled.bin(2), 2.5, "DynamicHistogram::bin failed on enabled histogram at 2");
    ASSERT_EQUAL(histogramEnabled.bin(3), 3.5, "DynamicHistogram::bin failed on enabled histogram at 3");

    std::vector<double> cdf = histogramEnabled.cdf();
    ASSERT_EQUAL(cdf.at(0), 0.25, "DynamicHistogram::cdf failed on enabled histogram at 0");
    ASSERT_EQUAL(cdf.at(1), 0.50, "DynamicHistogram::cdf failed on enabled histogram at 1");
    ASSERT_EQUAL(cdf.at(2), 0.75, "DynamicHistogram::cdf failed on enabled histogram at 2");
    ASSERT_EQUAL(cdf.at(3), 1.00, "DynamicHistogram::cdf failed on enabled histogram at 3");

    ASSERT(!histogramEnabled.contains(-1), "DynamicHistogram::contains failed on enabled histogram at -1");
    ASSERT(histogramEnabled.contains(0), "DynamicHistogram::contains failed on enabled histogram at 0");
    ASSERT(histogramEnabled.contains(4), "DynamicHistogram::contains failed on enabled histogram at 4");
    ASSERT(!histogramEnabled.contains(5), "DynamicHistogram::contains failed on enabled histogram at 5");

    ASSERT_EQUAL(histogramEnabled.count(), 4, "DynamicHistogram::count failed on enabled histogram");

    ASSERT(histogramEnabled.enabled(), "DynamicHistogram::enabled failed on enabled histogram");
    
    ASSERT_EQUAL(histogramEnabled.frequency(0), 2, "DynamicHistogram::frequency failed on enabled histogram at 0");
    ASSERT_EQUAL(histogramEnabled.frequency(1), 2, "DynamicHistogram::frequency failed on enabled histogram at 1");
    ASSERT_EQUAL(histogramEnabled.frequency(2), 2, "DynamicHistogram::frequency failed on enabled histogram at 2");
    ASSERT_EQUAL(histogramEnabled.frequency(3), 2, "DynamicHistogram::frequency failed on enabled histogram at 3");

    std::vector<double> frequencies = histogramEnabled.frequencies();
    ASSERT_EQUAL(frequencies.at(0), 2, "DynamicHistogram::frequencies failed on enabled histogram at 0");
    ASSERT_EQUAL(frequencies.at(1), 2, "DynamicHistogram::frequencies failed on enabled histogram at 1");
    ASSERT_EQUAL(frequencies.at(2), 2, "DynamicHistogram::frequencies failed on enabled histogram at 2");
    ASSERT_EQUAL(frequencies.at(3), 2, "DynamicHistogram::frequencies failed on enabled histogram at 3");

    ASSERT_EQUAL(histogramEnabled.index(0), 0, "DynamicHistogram::index failed on enabled histogram at 0");
    ASSERT_EQUAL(histogramEnabled.index(1), 1, "DynamicHistogram::index failed on enabled histogram at 1");
    ASSERT_EQUAL(histogramEnabled.index(2), 2, "DynamicHistogram::index failed on enabled histogram at 2");
    ASSERT_EQUAL(histogramEnabled.index(4), 3, "DynamicHistogram::index failed on enabled histogram at 4");

    ASSERT(histogramEnabled.initialized(), "DynamicHistogram::initialized failed on enabled histogram");
    ASSERT_EQUAL(histogramEnabled.max(), 4, "DynamicHistogram::max failed on enabled histogram");
    ASSERT_EQUAL(histogramEnabled.merges(), 0, "DynamicHistogram::merges failed on enabled histogram");
    ASSERT_EQUAL(histogramEnabled.min(), 0, "DynamicHistogram::min failed on enabled histogram");

    std::vector<double> pdf = histogramEnabled.pdf();
    ASSERT_EQUAL(pdf.at(0), 0.25, "DynamicHistogram::pdf failed on enabled histogram at 0");
    ASSERT_EQUAL(pdf.at(1), 0.25, "DynamicHistogram::pdf failed on enabled histogram at 1");
    ASSERT_EQUAL(pdf.at(2), 0.25, "DynamicHistogram::pdf failed on enabled histogram at 2");
    ASSERT_EQUAL(pdf.at(3), 0.25, "DynamicHistogram::pdf failed on enabled histogram at 3");

    ASSERT_EQUAL(histogramEnabled.total(), 8, "DynamicHistogram::total failed on enabled histogram");

    ASSERT_EQUAL(histogramEnabled.width(0), 1, "DynamicHistogram::width failed on enabled histogram at 0");
    ASSERT_EQUAL(histogramEnabled.width(1), 1, "DynamicHistogram::width failed on enabled histogram at 1");
    ASSERT_EQUAL(histogramEnabled.width(2), 1, "DynamicHistogram::width failed on enabled histogram at 2");
    ASSERT_EQUAL(histogramEnabled.width(4), 1, "DynamicHistogram::width failed on enabled histogram at 4");

}

/**
 * Test DataVector
 */
void test_datavector()
{

    DataVector src;
    src.push_back(DataPoint(1,true));
    src.push_back(DataPoint(2,true));
    src.push_back(DataPoint(3,true));

    std::vector<bool> mask(3,true);
    DataVector dst;
    
    dst.copy(src,mask);
    ASSERT_EQUAL(dst.at(0).value, 1, "DataVector::copy value failed on enabled data point at 0");
    ASSERT(dst.at(0).active, "DataVector::copy active failed on enabled data point at 0");
    ASSERT_EQUAL(dst.at(1).value, 2, "DataVector::copy value failed on enabled data point at 1");
    ASSERT(dst.at(1).active, "DataVector::copy active failed on enabled data point at 1");
    ASSERT_EQUAL(dst.at(2).value, 3, "DataVector::copy value failed on enabled data point at 2");
    ASSERT(dst.at(2).active, "DataVector::copy active failed on enabled data point at 2");

    src.deactivate();
    dst.clear();
    dst.copy(src,mask);
    ASSERT_EXCEPTION(dst.at(0), "DataVector::copy failed on disabled data point at 0");
    ASSERT_EXCEPTION(dst.at(1), "DataVector::copy failed on disabled data point at 1");
    ASSERT_EXCEPTION(dst.at(2), "DataVector::copy failed on disabled data point at 2");

}

/**
 * Test DataFilters
 */
void test_datafilters()
{

    std::vector<std::string> datas;
    datas.push_back("1");
    datas.push_back("2");
    datas.push_back("3");
    
    DataVector datav;
    datav.push_back(DataPoint(1,true));
    datav.push_back(DataPoint(2,true));
    datav.push_back(DataPoint(3,true));

    DataFilters filters;
    ASSERT(filters.isFiltered(datas), "DataFilters::isFiltered failed on string vector with no filter");
    ASSERT(filters.isFiltered(datav), "DataFilters::isFiltered failed on data vector with no filter");

    filters.addNumericFilter(0,1,1,true);
    ASSERT(filters.isFiltered(datav), "DataFilters::isFiltered failed on data vector with numeric accepting filter on 0");
    
    filters.addNumericFilter(1,4,4,false);
    ASSERT(filters.isFiltered(datav), "DataFilters::isFiltered failed on data vector with numeric rejection filter on 1");

    filters.addStringFilter(0,"1",true,true,true);
    ASSERT(filters.isFiltered(datas), "DataFilters::isFiltered failed on data vector with string accepting filter on 0");
    
    filters.addStringFilter(1,"4",true,true,false);
    ASSERT(filters.isFiltered(datas), "DataFilters::isFiltered failed on data vector with string rejection filter on 1");

}

/**
 * Test StatisticsTracker
 */
void test_statisticstracker()
{

    StatisticsTrackerOptions options;
    options.doCov = options.doMax = options.doMean = options.doMin = options.doVar = true;
    
    StatisticsTracker tracker("test", 0, 0, 0, 0, 0, options);
    tracker.update(1);
    tracker.update(2);
    tracker.update(3);
    tracker.update(4);
    tracker.update(5);

    ASSERT_EQUAL(tracker.getCount(), 5, "StatisticsTracker::getCount failed");
    ASSERT_EQUAL(tracker.getMinimum(), 1, "StatisticsTracker::getMin failed");
    ASSERT_EQUAL(tracker.getMean(), 3, "StatisticsTracker::getMean failed");
    ASSERT_EQUAL(tracker.getMaximum(), 5, "StatisticsTracker::getMax failed");
    ASSERT_EQUAL(tracker.getName(), "test", "StatisticsTracker::getName failed");
    ASSERT_EQUAL(tracker.getVariance(), 2, "StatisticsTracker::getVariance failed");
    
    tracker.setName("test2");
    ASSERT_EQUAL(tracker.getName(), "test2", "StatisticsTracker::setName failed");

}

/**
 * Test MultivariateTracker
 */
void test_multivariatetracker()
{

    DataVector data1, data2, data3;
    data1.push_back(DataPoint(1,true));
    data1.push_back(DataPoint(3,true));
    data2.push_back(DataPoint(2,true));
    data2.push_back(DataPoint(2,true));
    data3.push_back(DataPoint(3,true));
    data3.push_back(DataPoint(1,true));

    StatisticsTrackerOptions options;
    options.doCov = options.doMax = options.doMean = options.doMin = options.doVar = true;
    MultivariateTracker tracker(2, options);

    tracker.update(data1);
    tracker.update(data2);
    tracker.update(data3);

    double cov = 2.0/3.0;
    ASSERT_EQUAL(tracker.getCovariance(0,0), cov, "StatisticsTracker::getCovariance failed on (0,0)");
    ASSERT_EQUAL(tracker.getCovariance(0,1), -cov, "StatisticsTracker::getCovariance failed on (0,1)");
    ASSERT_EQUAL(tracker.getCovariance(1,0), -cov, "StatisticsTracker::getCovariance failed on (1,0)");
    ASSERT_EQUAL(tracker.getCovariance(1,1), cov, "StatisticsTracker::getCovariance failed on (1,1)");

}

/**
 * Test RowSampler uniform sampling
 */
void test_rowsampler_uniform()
{

    SamplerOptions options;
    options.mode = SamplerOptions::UNIFORM;
    options.step = 10;
    
    RowSampler sampler(options);
    
    int n = 100;
    for (int i = 0; i < n; i++)
    {
        bool isSampledTest = sampler.sample();
        bool isSampledTrue = (i % (options.step)) == 0;
        ASSERT_EQUAL(isSampledTest, isSampledTrue, "RowSampler::sample with uniform failed on " << i);
    }

}

/**
 * Test RowSampler random sampling
 */
void test_rowsampler_random()
{

    SamplerOptions options;
    options.mode = SamplerOptions::RANDOM;
    options.step = 10;
    
    RowSampler sampler(options);
    
    int n = 100;
    int numSamplesTest = 0;
    for (int i = 0; i < n; i++)
    {
        if (sampler.sample()) numSamplesTest++;
    }
    
    int numSamplesTrue = n/options.step;
    ASSERT_EQUAL(numSamplesTest, numSamplesTrue, "RowSampler::sample with random failed");

}

/**
 * Unit tester for clistats
 */
int
main(int argc,
     char * argv[])
{
    
    std::vector<void (*)()> tests;
    std::vector<std::string> errors;
    
    // ======================================================================
    // Define tests to run
    // ----------------------------------------------------------------------
    tests.push_back(test_stringparser_parsenumer);
    tests.push_back(test_stringparser_parsestatistic);
    tests.push_back(test_stringparser_replacen);
    tests.push_back(test_stringparser_tovalue);
    tests.push_back(test_stringparser_trim);
    tests.push_back(test_stringparser_updatewidthnumber);
    tests.push_back(test_stringparser_updatewidthstring);
    tests.push_back(test_stringsplitter);
    tests.push_back(test_stringsplitter_tointegers);
    tests.push_back(test_fixedsizecache);
    tests.push_back(test_dynamichistogram_disabled);
    tests.push_back(test_dynamichistogram_enabled);
    tests.push_back(test_datavector);
    tests.push_back(test_datafilters);
    tests.push_back(test_statisticstracker);
    tests.push_back(test_multivariatetracker);
    tests.push_back(test_rowsampler_uniform);
    tests.push_back(test_rowsampler_random);

    // ======================================================================
    // Define test metrics
    // ----------------------------------------------------------------------
    unsigned int numTestsTotal = tests.size();
    unsigned int numTestsPass = 0;
    unsigned int numTestsFail = 0;

    // ======================================================================
    // Display parameters
    // ----------------------------------------------------------------------
    unsigned int width = 70;
    std::string header(width,'=');
    std::string footer(width,'-');

    // ======================================================================
    // Display header
    // ----------------------------------------------------------------------
    std::cout << std::endl;
    std::cout << "        Testing " << ApplicationProperties::NAME  << " " << ApplicationProperties::VERSION_MAJOR << "." << ApplicationProperties::VERSION_MINOR << std::endl;
    std::cout << std::endl;

    std::cout << header << std::endl;
    std::cout << "    Running Tests (. = Pass, F = Pass)" << std::endl;
    std::cout << footer << std::endl;

    // ======================================================================
    // Run tests
    // ----------------------------------------------------------------------
    for (std::vector<void (*)()>::iterator test = tests.begin(); test != tests.end(); ++test)
    {
        std::string message = ".";
        try
        {
            (*test)();
        }
        catch (std::exception const & ex)
        {
            numTestsFail++;
            errors.push_back(ex.what());
            message = "F";
        }
        std::cout << message;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    
    numTestsPass = numTestsTotal - numTestsFail;

    // ======================================================================
    // Display test results
    // ----------------------------------------------------------------------
    std::cout << header << std::endl;
    std::cout << "    Testing results " << std::endl;
    std::cout << footer << std::endl;
    std::cout << "NUMBER TESTS TOTAL: " << numTestsTotal << std::endl;
    std::cout << "NUMBER TESTS PASS : " << numTestsPass << std::endl;
    std::cout << "NUMBER TESTS FAIL : " << numTestsFail << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    // Exit if no errors
    if (numTestsFail == 0)
    {
        std::cout << "OK" << std::endl;
        return 0;
    }

    // Print errors that occured
    std::cerr << header << std::endl;
    std::cerr << "    Testing Errors " << std::endl;
    std::cerr << footer << std::endl;
    unsigned int index = 0;
    for (std::vector<std::string>::iterator error = errors.begin(); error != errors.end(); ++error)
    {
        std::cerr << "Error #" << ++index << std::endl;
        std::cerr << *error << std::endl;
        std::cerr << footer << std::endl;
    }
    std::cerr << std::endl;

    return 1;
    
}
