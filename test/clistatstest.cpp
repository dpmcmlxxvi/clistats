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
            (void)(condition);                      \
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

    ASSERT_EQUAL(StringParser::parseStatistic<short>(1,1), "1", "StringParser::parseStatistic<short> failed on non-zero count");
    ASSERT_EQUAL(StringParser::parseStatistic<int>(1,1), "1", "StringParser::parseStatistic<int> failed on non-zero count");
    ASSERT_EQUAL(StringParser::parseStatistic<long>(1,1), "1", "StringParser::parseStatistic<long> failed on non-zero count");
    ASSERT_EQUAL(StringParser::parseStatistic<long long>(1,1), "1", "StringParser::parseStatistic<long long> failed on non-zero count");
    ASSERT_EQUAL(StringParser::parseStatistic<float>(1,(float)1.1), "1.100000", "StringParser::parseStatistic<float> failed on non-zero count");
    ASSERT_EQUAL(StringParser::parseStatistic<long double>(1,1.1), "1.100000", "StringParser::parseStatistic<long double> failed on non-zero count");

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

    int width = 0;
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

    int width = 0;
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

    const long count = 3;
    std::string valuesstring = "foo,bar,1";
    const char * valueschar[] = {"foo", "bar", "1"};

    StringSplitter splitter(valuesstring,",");
    for (long i = 0; i < count; i++)
    {
        ASSERT_EQUAL(splitter.at(i), valueschar[i], "StringSplitter::at failed on " << valueschar[i]);
    }

    std::vector<std::string> valuesGet(splitter.get());
    for (long i = 0; i < count; i++)
    {
        ASSERT_EQUAL(valueschar[i], valuesGet.at(i), "StringSplitter::get failed at " << valueschar[i]);
    }

    ASSERT_EQUAL(splitter.size(), count, "StringSplitter::size failed");

    std::vector<std::string> valuesTokens(splitter.tokens());
    for (long i = 0; i < count; i++)
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

    FixedSizeCache<int> cacheMinMax(3);
    cacheMinMax.add(3);
    cacheMinMax.add(1);
    cacheMinMax.add(2);
    ASSERT_EQUAL(cacheMinMax.max(), 3, "FixedSizeCache::max failed on non-empty cache");
    ASSERT_EQUAL(cacheMinMax.min(), 1, "FixedSizeCache::min failed on non-empty cache");

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
    ASSERT_EXCEPTION(histogramDisabled.width(), "DynamicHistogram::width failed on disabled histogram");

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

    ASSERT_EQUAL(histogramEnabled.width(), 1, "DynamicHistogram::width failed on enabled histogram at 0");

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

    DataFilters caseFilter;
    caseFilter.addStringFilter(0,"A",false,false,true);
    std::vector<std::string> matched;
    matched.push_back("abc");
    ASSERT(caseFilter.isFiltered(matched), "DataFilters::isFiltered failed on case-insensitive partial match");

    DataFilters numericMiss;
    numericMiss.addNumericFilter(0,100,200,true);
    ASSERT(!numericMiss.isFiltered(datav), "DataFilters::isFiltered failed to reject data vector with no matching numeric filter");

    DataFilters stringMiss;
    stringMiss.addStringFilter(0,"zzz",true,true,true);
    ASSERT(!stringMiss.isFiltered(datas), "DataFilters::isFiltered failed to reject string vector with no matching string filter");

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

    DataVector data4, data5, data6;
    data4.push_back(DataPoint(1,true));
    data4.push_back(DataPoint(2,true));
    data5.push_back(DataPoint(2,true));
    data5.push_back(DataPoint(4,true));
    data6.push_back(DataPoint(3,true));
    data6.push_back(DataPoint(6,true));

    MultivariateTracker slopeTracker(2, options);
    slopeTracker.update(data4);
    slopeTracker.update(data5);
    slopeTracker.update(data6);

    ASSERT_EQUAL(slopeTracker.getLeastSquaresSlope(0,1), 2.0, "MultivariateTracker::getLeastSquaresSlope failed");
    ASSERT_EQUAL(slopeTracker.getLeastSquaresOffset(0,1), 0.0, "MultivariateTracker::getLeastSquaresOffset failed");

    DataVector mismatched;
    mismatched.push_back(DataPoint(1,true));
    ASSERT(!slopeTracker.update(mismatched), "MultivariateTracker::update failed to reject data with mismatched dimensions");

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
 * Test DynamicHistogram with all identical values (degenerate range)
 * Robustness: must initialize and report counts without throw/empty histogram
 */
void test_dynamichistogram_degenerate()
{

    DynamicHistogramOptions options;
    options.enabled = true;
    options.binCount = 10;
    options.cacheSize = 8;

    DynamicHistogram<double> histogram(options);

    // Fill cache with identical values (triggers zero-range init on merge)
    for (int i = 0; i < options.cacheSize; i++)
    {
        histogram.add(0.0);
    }

    // Force merge / materialize histogram
    double total = histogram.total();
    ASSERT(total > 0, "DynamicHistogram degenerate: total should be > 0 after identical values");
    ASSERT_EQUAL(histogram.count(), options.binCount, "DynamicHistogram degenerate: bin count mismatch");

    // All mass should be present (sum of frequencies)
    std::vector<double> freq = histogram.frequencies();
    double sum = 0;
    for (std::vector<double>::size_type i = 0; i < freq.size(); i++)
    {
        sum += freq.at(i);
    }
    ASSERT_EQUAL(sum, (double) options.cacheSize, "DynamicHistogram degenerate: frequency sum should equal number of added points");

}

/**
 * Test DynamicHistogram expansion after degenerate init (zeros then a distant value)
 * Robustness: merge expansion must not throw on large range jump
 */
void test_dynamichistogram_expand_after_degenerate()
{

    DynamicHistogramOptions options;
    options.enabled = true;
    options.binCount = 10;
    options.cacheSize = 4;

    DynamicHistogram<double> histogram(options);

    for (int i = 0; i < options.cacheSize; i++)
    {
        histogram.add(0.0);
    }

    // Force initial merge on identical values
    ASSERT(histogram.total() > 0, "DynamicHistogram expand: expected positive total after zeros");

    // Distant value should expand histogram without throwing
    histogram.add(1.0);
    histogram.add(100.0);

    double total = histogram.total();
    ASSERT(total >= (double)(options.cacheSize + 2), "DynamicHistogram expand: total should include zeros and new points");

    // Accessors must remain usable
    std::vector<double> freq = histogram.frequencies();
    ASSERT_EQUAL((int) freq.size(), options.binCount, "DynamicHistogram expand: frequency vector size mismatch");

}

/**
 * Test MultivariateTracker covariance/correlation with inactive columns
 * Robustness: inactive entries must not corrupt covariance or throw
 */
void test_multivariatetracker_inactive()
{

    StatisticsTrackerOptions options;
    options.doCov = options.doMax = options.doMean = options.doMin = options.doVar = true;
    MultivariateTracker tracker(2, options);

    DataVector row1;
    row1.push_back(DataPoint(1.0, true));
    row1.push_back(DataPoint(2.0, false));  // inactive

    DataVector row2;
    row2.push_back(DataPoint(3.0, true));
    row2.push_back(DataPoint(4.0, true));

    DataVector row3;
    row3.push_back(DataPoint(5.0, true));
    row3.push_back(DataPoint(6.0, true));

    ASSERT(tracker.update(row1), "MultivariateTracker inactive: update row1 failed");
    ASSERT(tracker.update(row2), "MultivariateTracker inactive: update row2 failed");
    ASSERT(tracker.update(row3), "MultivariateTracker inactive: update row3 failed");

    // Column 0 saw 3 values; column 1 saw 2 (first inactive)
    ASSERT_EQUAL(tracker.getCount(0), 3, "MultivariateTracker inactive: count column 0");
    ASSERT_EQUAL(tracker.getCount(1), 2, "MultivariateTracker inactive: count column 1");

    // Covariance diagonal should be finite and non-negative
    double v0 = tracker.getCovariance(0, 0);
    double v1 = tracker.getCovariance(1, 1);
    ASSERT(v0 >= 0.0, "MultivariateTracker inactive: variance(0,0) should be >= 0");
    ASSERT(v1 >= 0.0, "MultivariateTracker inactive: variance(1,1) should be >= 0");

    // Correlation must not NaN when variances are positive
    double corr = tracker.getCorrelation(0, 1);
    ASSERT(corr == corr, "MultivariateTracker inactive: correlation should not be NaN");
    ASSERT(corr >= -1.0 && corr <= 1.0, "MultivariateTracker inactive: correlation out of [-1,1]");

}

/**
 * Test correlation when a column has zero variance
 * Robustness: must not return NaN
 */
void test_multivariatetracker_zero_variance_correlation()
{

    StatisticsTrackerOptions options;
    options.doCov = options.doMax = options.doMean = options.doMin = options.doVar = true;
    MultivariateTracker tracker(2, options);

    DataVector row1;
    row1.push_back(DataPoint(1.0, true));
    row1.push_back(DataPoint(7.0, true));

    DataVector row2;
    row2.push_back(DataPoint(2.0, true));
    row2.push_back(DataPoint(7.0, true));  // constant column

    DataVector row3;
    row3.push_back(DataPoint(3.0, true));
    row3.push_back(DataPoint(7.0, true));

    tracker.update(row1);
    tracker.update(row2);
    tracker.update(row3);

    double corr = tracker.getCorrelation(0, 1);
    ASSERT(corr == corr, "MultivariateTracker zero-variance: correlation should not be NaN");
    ASSERT_EQUAL(corr, 0.0, "MultivariateTracker zero-variance: correlation should be 0");

}

/**
 * Test StringSplitter::toIntegers range and single values used by column filters
 * Robustness: invalid ranges must fail cleanly; valid ranges must expand
 */
void test_stringsplitter_tointegers_filter_edges()
{

    std::vector<int> range;

    range.clear();
    ASSERT(StringSplitter::toIntegers("3:5", range), "toIntegers filter edges: valid range 3:5 should parse");
    ASSERT_EQUAL((int) range.size(), 3, "toIntegers filter edges: 3:5 size");
    ASSERT_EQUAL(range.at(0), 3, "toIntegers filter edges: 3:5 start");
    ASSERT_EQUAL(range.at(2), 5, "toIntegers filter edges: 3:5 end");

    range.clear();
    ASSERT(StringSplitter::toIntegers("2", range), "toIntegers filter edges: single value should parse");
    ASSERT_EQUAL((int) range.size(), 1, "toIntegers filter edges: single size");
    ASSERT_EQUAL(range.at(0), 2, "toIntegers filter edges: single value");

    range.clear();
    ASSERT(!StringSplitter::toIntegers("a:b", range), "toIntegers filter edges: non-numeric should fail");
    ASSERT(!StringSplitter::toIntegers("1:2:3", range), "toIntegers filter edges: too many tokens should fail");

    range.clear();
    ASSERT(StringSplitter::toIntegers("5:3", range), "toIntegers filter edges: reversed range should parse");
    ASSERT_EQUAL((int) range.size(), 3, "toIntegers filter edges: 5:3 size");
    ASSERT_EQUAL(range.at(0), 3, "toIntegers filter edges: 5:3 start");
    ASSERT_EQUAL(range.at(2), 5, "toIntegers filter edges: 5:3 end");

}

/**
 * Build a mutable argument array from a list of strings for CommandLineParser tests
 * @param[in] args Argument values, args.at(0) is treated as the executable name
 * @return Vector of char pointers into args suitable for use as argv
 */
std::vector<char *>
buildArgs(std::vector<std::string> & args)
{
    std::vector<char *> argv;
    for (std::vector<std::string>::iterator it = args.begin(); it != args.end(); ++it)
    {
        argv.push_back(&(*it)[0]);
    }
    return argv;
}

/**
 * Test CommandLineParser input/output/format flags
 */
void test_commandlineparser_flags_io()
{

    std::string fileInput = "test_commandlineparser_flags_io_input.csv";
    std::string fileOutput = "test_commandlineparser_flags_io_output.csv";
    std::ofstream input(fileInput.c_str());
    input << "1,2" << std::endl;
    input.close();

    std::vector<std::string> args;
    args.push_back("clistats");
    args.push_back("-i"); args.push_back(fileInput);
    args.push_back("-o"); args.push_back(fileOutput);
    args.push_back("-d"); args.push_back(";");
    args.push_back("-c"); args.push_back("#");
    args.push_back("-a"); args.push_back("44");
    args.push_back("-b");
    args.push_back("-s"); args.push_back("2");
    args.push_back("-k"); args.push_back("5");
    args.push_back("-r");
    args.push_back("-t"); args.push_back("1");
    std::vector<char *> argv = buildArgs(args);

    CommandLineParser parser((int) argv.size(), &argv[0]);

    ASSERT_EQUAL(parser.options.fileInput, fileInput, "CommandLineParser -i failed to set input file");
    ASSERT_EQUAL(parser.options.fileOutput, fileOutput, "CommandLineParser -o failed to set output file");
    ASSERT_EQUAL(parser.options.delimiter, ";,", "CommandLineParser -d/-a failed to append delimiter characters");
    ASSERT_EQUAL(parser.options.comment, "#", "CommandLineParser -c failed to set comment character");
    ASSERT(parser.options.blankEOF, "CommandLineParser -b failed to set blank EOF flag");
    ASSERT_EQUAL(parser.options.numLinesToSkip, (unsigned int) 2, "CommandLineParser -s failed to set number of lines to skip");
    ASSERT_EQUAL(parser.options.numLinesToKeep, (unsigned int) 5, "CommandLineParser -k failed to set number of lines to keep");
    ASSERT(parser.options.removeDuplicates, "CommandLineParser -r failed to set remove duplicates flag");
    ASSERT_EQUAL(parser.options.headerRow, (unsigned int) 1, "CommandLineParser -t failed to set header row");

    std::remove(fileInput.c_str());
    std::remove(fileOutput.c_str());

}

/**
 * Test CommandLineParser column and row filter flags
 */
void test_commandlineparser_flags_filters()
{

    std::vector<std::string> args;
    args.push_back("clistats");
    args.push_back("-fc"); args.push_back("1,3");
    args.push_back("-fn"); args.push_back("2,0,10");
    args.push_back("-fs"); args.push_back("4,abc");
    std::vector<char *> argv = buildArgs(args);

    CommandLineParser parser((int) argv.size(), &argv[0]);

    ASSERT(parser.options.filterColumns.at(0), "CommandLineParser -fc failed to enable column 1");
    ASSERT(!parser.options.filterColumns.at(1), "CommandLineParser -fc failed to leave column 2 disabled");
    ASSERT(parser.options.filterColumns.at(2), "CommandLineParser -fc failed to enable column 3");

}

/**
 * Test CommandLineParser sampling and parsing mode flags
 */
void test_commandlineparser_flags_sampling()
{

    std::vector<std::string> args;
    args.push_back("clistats");
    args.push_back("-rs"); args.push_back("2");
    args.push_back("-se"); args.push_back("7");
    args.push_back("-su"); args.push_back("3");
    args.push_back("-st");
    std::vector<char *> argv = buildArgs(args);

    CommandLineParser parser((int) argv.size(), &argv[0]);

    ASSERT_EQUAL(parser.options.numLinesToReshape, (unsigned int) 2, "CommandLineParser -rs failed to set reshape count");
    ASSERT_EQUAL(parser.options.seed, (unsigned int) 7, "CommandLineParser -se failed to set seed");
    ASSERT_EQUAL(parser.options.sampling.mode, SamplerOptions::UNIFORM, "CommandLineParser -su failed to set uniform sampling mode");
    ASSERT_EQUAL(parser.options.sampling.step, (unsigned int) 3, "CommandLineParser -su failed to set sampling step");
    ASSERT(parser.options.strictParsing, "CommandLineParser -st failed to set strict parsing flag");

    std::vector<std::string> argsRandom;
    argsRandom.push_back("clistats");
    argsRandom.push_back("-sr"); argsRandom.push_back("5");
    std::vector<char *> argvRandom = buildArgs(argsRandom);

    CommandLineParser parserRandom((int) argvRandom.size(), &argvRandom[0]);
    ASSERT_EQUAL(parserRandom.options.sampling.mode, SamplerOptions::RANDOM, "CommandLineParser -sr failed to set random sampling mode");
    ASSERT_EQUAL(parserRandom.options.sampling.step, (unsigned int) 5, "CommandLineParser -sr failed to set sampling step");

}

/**
 * Test CommandLineParser output selection flags
 */
void test_commandlineparser_flags_output()
{

    std::vector<std::string> args;
    args.push_back("clistats");
    args.push_back("-v"); args.push_back("2");
    args.push_back("-cv");
    args.push_back("-cr");
    args.push_back("-fd");
    args.push_back("-hg"); args.push_back("5,200");
    args.push_back("-lo");
    args.push_back("-ls");
    args.push_back("-ss");
    std::vector<char *> argv = buildArgs(args);

    CommandLineParser parser((int) argv.size(), &argv[0]);

    ASSERT_EQUAL((int) parser.options.verboseLevel, 2, "CommandLineParser -v failed to set verbose level");
    ASSERT(parser.options.showCovariance, "CommandLineParser -cv failed to set covariance flag");
    ASSERT(parser.options.showCorrelation, "CommandLineParser -cr failed to set correlation flag");
    ASSERT(parser.options.showFilteredData, "CommandLineParser -fd failed to set filtered data flag");
    ASSERT(parser.options.showHistogram, "CommandLineParser -hg failed to set histogram flag");
    ASSERT_EQUAL(parser.options.statisticsOptions.histogramOptions.binCount, 5, "CommandLineParser -hg failed to set bin count");
    ASSERT_EQUAL(parser.options.statisticsOptions.histogramOptions.cacheSize, 200, "CommandLineParser -hg failed to set cache size");
    ASSERT(parser.options.showLeastSquaresOffset, "CommandLineParser -lo failed to set least squares offset flag");
    ASSERT(parser.options.showLeastSquaresSlope, "CommandLineParser -ls failed to set least squares slope flag");
    ASSERT(parser.options.showStatistics, "CommandLineParser -ss failed to set statistics flag");

}

/**
 * Test CommandLineParser help and version flags
 */
void test_commandlineparser_help_version()
{

    std::vector<std::string> argsHelp;
    argsHelp.push_back("clistats");
    argsHelp.push_back("-h");
    std::vector<char *> argvHelp = buildArgs(argsHelp);

    CommandLineParser parserHelp((int) argvHelp.size(), &argvHelp[0]);
    ASSERT(parserHelp.showUsage(), "CommandLineParser -h failed to set usage flag");

    std::vector<std::string> argsVersion;
    argsVersion.push_back("clistats");
    argsVersion.push_back("-V");
    std::vector<char *> argvVersion = buildArgs(argsVersion);

    CommandLineParser parserVersion((int) argvVersion.size(), &argvVersion[0]);
    ASSERT(parserVersion.showVersion(), "CommandLineParser -V failed to set version flag");

    std::ostringstream usage;
    std::streambuf * oldUsageBuf = std::cout.rdbuf(usage.rdbuf());
    CommandLineParser::printUsage();
    std::cout.rdbuf(oldUsageBuf);
    ASSERT(usage.str().find("SYNOPSIS") != std::string::npos, "CommandLineParser::printUsage failed to print usage text");

    std::ostringstream version;
    std::streambuf * oldVersionBuf = std::cout.rdbuf(version.rdbuf());
    CommandLineParser::printVersion();
    std::cout.rdbuf(oldVersionBuf);
    ASSERT_EQUAL(version.str(), ApplicationProperties::version() + "\n", "CommandLineParser::printVersion failed to print version text");

}

/**
 * Test CommandLineParser rejects unrecognized and invalid flag values
 */
void test_commandlineparser_invalid()
{

    std::vector<std::string> argsUnknown;
    argsUnknown.push_back("clistats");
    argsUnknown.push_back("--bogus");
    std::vector<char *> argvUnknown = buildArgs(argsUnknown);
    ASSERT_EXCEPTION(CommandLineParser(((int) argvUnknown.size()), &argvUnknown[0]), "CommandLineParser failed to reject unrecognized flag");

    std::vector<std::string> argsBadVerbose;
    argsBadVerbose.push_back("clistats");
    argsBadVerbose.push_back("-v"); argsBadVerbose.push_back("99");
    std::vector<char *> argvBadVerbose = buildArgs(argsBadVerbose);
    ASSERT_EXCEPTION(CommandLineParser((int) argvBadVerbose.size(), &argvBadVerbose[0]), "CommandLineParser failed to reject invalid verbose level");

    std::vector<std::string> argsBadInput;
    argsBadInput.push_back("clistats");
    argsBadInput.push_back("-i"); argsBadInput.push_back("does_not_exist.csv");
    std::vector<char *> argvBadInput = buildArgs(argsBadInput);
    ASSERT_EXCEPTION(CommandLineParser((int) argvBadInput.size(), &argvBadInput[0]), "CommandLineParser failed to reject missing input file");

    std::vector<std::string> argsSameDelimiterComment;
    argsSameDelimiterComment.push_back("clistats");
    argsSameDelimiterComment.push_back("-d"); argsSameDelimiterComment.push_back("#");
    argsSameDelimiterComment.push_back("-c"); argsSameDelimiterComment.push_back("#");
    std::vector<char *> argvSameDelimiterComment = buildArgs(argsSameDelimiterComment);
    ASSERT_EXCEPTION(CommandLineParser((int) argvSameDelimiterComment.size(), &argvSameDelimiterComment[0]), "CommandLineParser failed to reject matching delimiter and comment");

    std::vector<std::string> argsOutput;
    argsOutput.push_back("clistats");
    argsOutput.push_back("-o"); argsOutput.push_back("/no_such_directory_xyz/out.csv");
    std::vector<char *> argvOutput = buildArgs(argsOutput);
    ASSERT_EXCEPTION(CommandLineParser((int) argvOutput.size(), &argvOutput[0]), "CommandLineParser failed to reject unwritable output file");

    std::vector<std::string> argsAscii;
    argsAscii.push_back("clistats");
    argsAscii.push_back("-a"); argsAscii.push_back("200");
    std::vector<char *> argvAscii = buildArgs(argsAscii);
    ASSERT_EXCEPTION(CommandLineParser((int) argvAscii.size(), &argvAscii[0]), "CommandLineParser failed to reject out-of-range ASCII code");

    std::vector<std::string> argsComment;
    argsComment.push_back("clistats");
    argsComment.push_back("-c"); argsComment.push_back("##");
    std::vector<char *> argvComment = buildArgs(argsComment);
    ASSERT_EXCEPTION(CommandLineParser((int) argvComment.size(), &argvComment[0]), "CommandLineParser failed to reject multi-character comment");

    std::vector<std::string> argsSkip;
    argsSkip.push_back("clistats");
    argsSkip.push_back("-s"); argsSkip.push_back("-1");
    std::vector<char *> argvSkip = buildArgs(argsSkip);
    ASSERT_EXCEPTION(CommandLineParser((int) argvSkip.size(), &argvSkip[0]), "CommandLineParser failed to reject negative skip count");

    std::vector<std::string> argsKeep;
    argsKeep.push_back("clistats");
    argsKeep.push_back("-k"); argsKeep.push_back("-1");
    std::vector<char *> argvKeep = buildArgs(argsKeep);
    ASSERT_EXCEPTION(CommandLineParser((int) argvKeep.size(), &argvKeep[0]), "CommandLineParser failed to reject negative keep count");

    std::vector<std::string> argsTitles;
    argsTitles.push_back("clistats");
    argsTitles.push_back("-t"); argsTitles.push_back("0");
    std::vector<char *> argvTitles = buildArgs(argsTitles);
    ASSERT_EXCEPTION(CommandLineParser((int) argvTitles.size(), &argvTitles[0]), "CommandLineParser failed to reject non-positive header row");

    std::vector<std::string> argsFcBadValue;
    argsFcBadValue.push_back("clistats");
    argsFcBadValue.push_back("-fc"); argsFcBadValue.push_back("a");
    std::vector<char *> argvFcBadValue = buildArgs(argsFcBadValue);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFcBadValue.size(), &argvFcBadValue[0]), "CommandLineParser failed to reject non-numeric filter column");

    std::vector<std::string> argsFcZero;
    argsFcZero.push_back("clistats");
    argsFcZero.push_back("-fc"); argsFcZero.push_back("0");
    std::vector<char *> argvFcZero = buildArgs(argsFcZero);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFcZero.size(), &argvFcZero[0]), "CommandLineParser failed to reject non-positive filter column");

    std::vector<std::string> argsFnEmpty;
    argsFnEmpty.push_back("clistats");
    argsFnEmpty.push_back("-fn");
    std::vector<char *> argvFnEmpty = buildArgs(argsFnEmpty);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFnEmpty.size(), &argvFnEmpty[0]), "CommandLineParser failed to reject empty numeric filter");

    std::vector<std::string> argsFnFormat;
    argsFnFormat.push_back("clistats");
    argsFnFormat.push_back("-fn"); argsFnFormat.push_back("1,2");
    std::vector<char *> argvFnFormat = buildArgs(argsFnFormat);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFnFormat.size(), &argvFnFormat[0]), "CommandLineParser failed to reject numeric filter with too few tokens");

    std::vector<std::string> argsFnColumn;
    argsFnColumn.push_back("clistats");
    argsFnColumn.push_back("-fn"); argsFnColumn.push_back("a,0,10");
    std::vector<char *> argvFnColumn = buildArgs(argsFnColumn);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFnColumn.size(), &argvFnColumn[0]), "CommandLineParser failed to reject non-numeric numeric filter column");

    std::vector<std::string> argsFnMin;
    argsFnMin.push_back("clistats");
    argsFnMin.push_back("-fn"); argsFnMin.push_back("1,abc,10");
    std::vector<char *> argvFnMin = buildArgs(argsFnMin);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFnMin.size(), &argvFnMin[0]), "CommandLineParser failed to reject non-numeric numeric filter minimum");

    std::vector<std::string> argsFnMax;
    argsFnMax.push_back("clistats");
    argsFnMax.push_back("-fn"); argsFnMax.push_back("1,0,abc");
    std::vector<char *> argvFnMax = buildArgs(argsFnMax);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFnMax.size(), &argvFnMax[0]), "CommandLineParser failed to reject non-numeric numeric filter maximum");

    std::vector<std::string> argsFnAccept;
    argsFnAccept.push_back("clistats");
    argsFnAccept.push_back("-fn"); argsFnAccept.push_back("1,0,10,x");
    std::vector<char *> argvFnAccept = buildArgs(argsFnAccept);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFnAccept.size(), &argvFnAccept[0]), "CommandLineParser failed to reject invalid numeric filter accept option");

    std::vector<std::string> argsFnInf;
    argsFnInf.push_back("clistats");
    argsFnInf.push_back("-fn"); argsFnInf.push_back("1,-inf,inf,a");
    std::vector<char *> argvFnInf = buildArgs(argsFnInf);
    CommandLineParser parserFnInf((int) argvFnInf.size(), &argvFnInf[0]);
    DataVector infData;
    infData.push_back(DataPoint(std::numeric_limits<double>::max(),true));
    ASSERT(parserFnInf.options.filterRows.isFiltered(infData), "CommandLineParser -fn failed to parse -inf/inf bounds and explicit accept option");

    std::vector<std::string> argsFnInfReversed;
    argsFnInfReversed.push_back("clistats");
    argsFnInfReversed.push_back("-fn"); argsFnInfReversed.push_back("1,inf,-inf");
    std::vector<char *> argvFnInfReversed = buildArgs(argsFnInfReversed);
    CommandLineParser parserFnInfReversed((int) argvFnInfReversed.size(), &argvFnInfReversed[0]);
    ASSERT(!parserFnInfReversed.options.filterRows.isFiltered(infData), "CommandLineParser -fn failed to parse inf minimum and -inf maximum bounds");

    std::vector<std::string> argsFsFormat;
    argsFsFormat.push_back("clistats");
    argsFsFormat.push_back("-fs"); argsFsFormat.push_back("1");
    std::vector<char *> argvFsFormat = buildArgs(argsFsFormat);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFsFormat.size(), &argvFsFormat[0]), "CommandLineParser failed to reject string filter with too few tokens");

    std::vector<std::string> argsFsColumn;
    argsFsColumn.push_back("clistats");
    argsFsColumn.push_back("-fs"); argsFsColumn.push_back("a,pattern");
    std::vector<char *> argvFsColumn = buildArgs(argsFsColumn);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFsColumn.size(), &argvFsColumn[0]), "CommandLineParser failed to reject non-numeric string filter column");

    std::vector<std::string> argsFsMatch;
    argsFsMatch.push_back("clistats");
    argsFsMatch.push_back("-fs"); argsFsMatch.push_back("1,pattern,x");
    std::vector<char *> argvFsMatch = buildArgs(argsFsMatch);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFsMatch.size(), &argvFsMatch[0]), "CommandLineParser failed to reject invalid string filter matching option");

    std::vector<std::string> argsFsSensitivity;
    argsFsSensitivity.push_back("clistats");
    argsFsSensitivity.push_back("-fs"); argsFsSensitivity.push_back("1,pattern,e,x");
    std::vector<char *> argvFsSensitivity = buildArgs(argsFsSensitivity);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFsSensitivity.size(), &argvFsSensitivity[0]), "CommandLineParser failed to reject invalid string filter case sensitivity option");

    std::vector<std::string> argsFsAccept;
    argsFsAccept.push_back("clistats");
    argsFsAccept.push_back("-fs"); argsFsAccept.push_back("1,pattern,e,s,x");
    std::vector<char *> argvFsAccept = buildArgs(argsFsAccept);
    ASSERT_EXCEPTION(CommandLineParser((int) argvFsAccept.size(), &argvFsAccept[0]), "CommandLineParser failed to reject invalid string filter accept option");

    std::vector<std::string> argsFsExplicitAccept;
    argsFsExplicitAccept.push_back("clistats");
    argsFsExplicitAccept.push_back("-fs"); argsFsExplicitAccept.push_back("1,pattern,e,s,a");
    std::vector<char *> argvFsExplicitAccept = buildArgs(argsFsExplicitAccept);
    CommandLineParser parserFsExplicitAccept((int) argvFsExplicitAccept.size(), &argvFsExplicitAccept[0]);
    std::vector<std::string> patternData;
    patternData.push_back("pattern");
    ASSERT(parserFsExplicitAccept.options.filterRows.isFiltered(patternData), "CommandLineParser -fs failed to parse explicit accept option");

    std::vector<std::string> argsReshape;
    argsReshape.push_back("clistats");
    argsReshape.push_back("-rs"); argsReshape.push_back("0");
    std::vector<char *> argvReshape = buildArgs(argsReshape);
    ASSERT_EXCEPTION(CommandLineParser((int) argvReshape.size(), &argvReshape[0]), "CommandLineParser failed to reject non-positive reshape count");

    std::vector<std::string> argsSeed;
    argsSeed.push_back("clistats");
    argsSeed.push_back("-se"); argsSeed.push_back("-1");
    std::vector<char *> argvSeed = buildArgs(argsSeed);
    ASSERT_EXCEPTION(CommandLineParser((int) argvSeed.size(), &argvSeed[0]), "CommandLineParser failed to reject negative seed value");

    std::vector<std::string> argsRandomStep;
    argsRandomStep.push_back("clistats");
    argsRandomStep.push_back("-sr"); argsRandomStep.push_back("0");
    std::vector<char *> argvRandomStep = buildArgs(argsRandomStep);
    ASSERT_EXCEPTION(CommandLineParser((int) argvRandomStep.size(), &argvRandomStep[0]), "CommandLineParser failed to reject non-positive random sampling step");

    std::vector<std::string> argsUniformStep;
    argsUniformStep.push_back("clistats");
    argsUniformStep.push_back("-su"); argsUniformStep.push_back("0");
    std::vector<char *> argvUniformStep = buildArgs(argsUniformStep);
    ASSERT_EXCEPTION(CommandLineParser((int) argvUniformStep.size(), &argvUniformStep[0]), "CommandLineParser failed to reject non-positive uniform sampling step");

    std::vector<std::string> argsMultipleA;
    argsMultipleA.push_back("clistats");
    argsMultipleA.push_back("-su"); argsMultipleA.push_back("2");
    argsMultipleA.push_back("-sr"); argsMultipleA.push_back("3");
    std::vector<char *> argvMultipleA = buildArgs(argsMultipleA);
    ASSERT_EXCEPTION(CommandLineParser((int) argvMultipleA.size(), &argvMultipleA[0]), "CommandLineParser failed to reject random sampling after uniform sampling already set");

    std::vector<std::string> argsMultipleB;
    argsMultipleB.push_back("clistats");
    argsMultipleB.push_back("-sr"); argsMultipleB.push_back("2");
    argsMultipleB.push_back("-su"); argsMultipleB.push_back("3");
    std::vector<char *> argvMultipleB = buildArgs(argsMultipleB);
    ASSERT_EXCEPTION(CommandLineParser((int) argvMultipleB.size(), &argvMultipleB[0]), "CommandLineParser failed to reject uniform sampling after random sampling already set");

    std::vector<std::string> argsHgMissing;
    argsHgMissing.push_back("clistats");
    argsHgMissing.push_back("-hg");
    std::vector<char *> argvHgMissing = buildArgs(argsHgMissing);
    ASSERT_EXCEPTION(CommandLineParser((int) argvHgMissing.size(), &argvHgMissing[0]), "CommandLineParser failed to reject missing histogram options");

    std::vector<std::string> argsHgExtra;
    argsHgExtra.push_back("clistats");
    argsHgExtra.push_back("-hg"); argsHgExtra.push_back("1,2,3");
    std::vector<char *> argvHgExtra = buildArgs(argsHgExtra);
    ASSERT_EXCEPTION(CommandLineParser((int) argvHgExtra.size(), &argvHgExtra[0]), "CommandLineParser failed to reject histogram options with too many tokens");

    std::vector<std::string> argsHgBins;
    argsHgBins.push_back("clistats");
    argsHgBins.push_back("-hg"); argsHgBins.push_back("0,200");
    std::vector<char *> argvHgBins = buildArgs(argsHgBins);
    ASSERT_EXCEPTION(CommandLineParser((int) argvHgBins.size(), &argvHgBins[0]), "CommandLineParser failed to reject non-positive bin count");

    std::vector<std::string> argsHgCache;
    argsHgCache.push_back("clistats");
    argsHgCache.push_back("-hg"); argsHgCache.push_back("5,50");
    std::vector<char *> argvHgCache = buildArgs(argsHgCache);
    ASSERT_EXCEPTION(CommandLineParser((int) argvHgCache.size(), &argvHgCache[0]), "CommandLineParser failed to reject undersized histogram cache");

}

/**
 * Test program driver flag parsing, run, and display paths
 */
void test_pdrv()
{

    std::string fileInput = "test_pdrv_input.csv";
    std::ofstream input(fileInput.c_str());
    input << "1,2" << std::endl;
    input << "3,4" << std::endl;
    input.close();

    std::string fileOutput = "test_pdrv_output.csv";
    std::vector<std::string> argsSuccess;
    argsSuccess.push_back("clistats");
    argsSuccess.push_back("-i"); argsSuccess.push_back(fileInput);
    argsSuccess.push_back("-o"); argsSuccess.push_back(fileOutput);
    argsSuccess.push_back("-d"); argsSuccess.push_back(",");
    std::vector<char *> argvSuccess = buildArgs(argsSuccess);
    ASSERT_EQUAL(pdrv((int) argvSuccess.size(), &argvSuccess[0]), (int) AppStatus::SUCCESS, "pdrv failed on valid input");
    std::remove(fileOutput.c_str());

    std::vector<std::string> argsHelp;
    argsHelp.push_back("clistats");
    argsHelp.push_back("-h");
    std::vector<char *> argvHelp = buildArgs(argsHelp);
    std::ostringstream usage;
    std::streambuf * oldUsageBuf = std::cout.rdbuf(usage.rdbuf());
    int statusHelp = pdrv((int) argvHelp.size(), &argvHelp[0]);
    std::cout.rdbuf(oldUsageBuf);
    ASSERT_EQUAL(statusHelp, (int) AppStatus::SUCCESS, "pdrv failed on -h flag");
    ASSERT(usage.str().find("SYNOPSIS") != std::string::npos, "pdrv failed to print usage on -h flag");

    std::vector<std::string> argsVersion;
    argsVersion.push_back("clistats");
    argsVersion.push_back("-V");
    std::vector<char *> argvVersion = buildArgs(argsVersion);
    std::ostringstream version;
    std::streambuf * oldVersionBuf = std::cout.rdbuf(version.rdbuf());
    int statusVersion = pdrv((int) argvVersion.size(), &argvVersion[0]);
    std::cout.rdbuf(oldVersionBuf);
    ASSERT_EQUAL(statusVersion, (int) AppStatus::SUCCESS, "pdrv failed on -V flag");
    ASSERT_EQUAL(version.str(), ApplicationProperties::version() + "\n", "pdrv failed to print version on -V flag");

    std::vector<std::string> argsBad;
    argsBad.push_back("clistats");
    argsBad.push_back("--bogus");
    std::vector<char *> argvBad = buildArgs(argsBad);
    std::ostringstream fatal;
    std::streambuf * oldCerrBuf = std::cerr.rdbuf(fatal.rdbuf());
    int statusBad = pdrv((int) argvBad.size(), &argvBad[0]);
    std::cerr.rdbuf(oldCerrBuf);
    ASSERT_EQUAL(statusBad, (int) AppStatus::FAILED_PARSING, "pdrv failed to report parsing failure on unrecognized flag");

    std::remove(fileInput.c_str());

}

/**
 * Test DynamicHistogram order statistic interpolation
 */
void test_dynamichistogram_order()
{

    DynamicHistogramOptions options;
    options.enabled = true;
    options.binCount = 4;
    options.cacheSize = 4;

    DynamicHistogram<double> histogram(options);
    histogram.add(0.0);
    histogram.add(1.0);
    histogram.add(2.0);
    histogram.add(4.0);

    ASSERT_EQUAL(histogram.order(0), histogram.bin(0), "DynamicHistogram::order failed on lower probability bound");
    ASSERT_EQUAL(histogram.order((int) histogram.total()), histogram.bin(3), "DynamicHistogram::order failed on upper probability bound");

    double interpolated = histogram.order(2);
    ASSERT(interpolated >= histogram.bin(0) && interpolated <= histogram.bin(3), "DynamicHistogram::order failed to interpolate within histogram range");

    DynamicHistogramOptions optionsFlat;
    optionsFlat.enabled = true;
    optionsFlat.binCount = 4;
    optionsFlat.cacheSize = 4;

    DynamicHistogram<double> histogramFlat(optionsFlat);
    histogramFlat.add(0.0);
    histogramFlat.add(0.0);
    histogramFlat.add(1.0);
    histogramFlat.add(10.0);

    double midpoint = histogramFlat.bin(0) + 0.5 * (histogramFlat.bin(1) - histogramFlat.bin(0));
    ASSERT_EQUAL(histogramFlat.order(3), midpoint, "DynamicHistogram::order failed to interpolate flat CDF region");

    DynamicHistogram<double> * histogramHeap = new DynamicHistogram<double>(options);
    histogramHeap->add(1.0);
    delete histogramHeap;

}

/**
 * Test StatisticsWriter matrix, histogram, and data output
 */
void test_statisticswriter_output()
{

    DataVector data1, data2;
    data1.push_back(DataPoint(1,true));
    data1.push_back(DataPoint(3,true));
    data2.push_back(DataPoint(2,true));
    data2.push_back(DataPoint(2,true));

    StatisticsTrackerOptions options;
    options.doCov = options.doMax = options.doMean = options.doMin = options.doVar = true;
    options.histogramOptions.enabled = true;
    options.histogramOptions.binCount = 2;
    options.histogramOptions.cacheSize = 2;
    MultivariateTracker tracker(2, options);
    tracker.update(data1);
    tracker.update(data2);

    StatisticsWriter writer;

    std::ostringstream matrix;
    writer.writeMatrix<double>(matrix, tracker, "Covariance", &MultivariateTracker::getCovariance);
    ASSERT(matrix.str().find("Covariance") != std::string::npos, "StatisticsWriter::writeMatrix failed to write title");

    std::ostringstream histograms;
    writer.writeHistograms(histograms, tracker);
    ASSERT(histograms.str().find("Histogram") != std::string::npos, "StatisticsWriter::writeHistograms failed to write histogram section");

    std::ostringstream emptyHistograms;
    MultivariateTracker emptyTracker(0, options);
    writer.writeHistograms(emptyHistograms, emptyTracker);
    ASSERT(emptyHistograms.str().find("No valid data") != std::string::npos, "StatisticsWriter::writeHistograms failed to report no data on empty tracker");

    std::ostringstream matrixEmpty;
    writer.writeMatrix<double>(matrixEmpty, emptyTracker, "Covariance", &MultivariateTracker::getCovariance);
    ASSERT(matrixEmpty.str().find("No Data") != std::string::npos, "StatisticsWriter::writeMatrix failed to report no data on empty tracker");

    DataVector data;
    data.push_back(DataPoint(1,true));
    data.push_back(DataPoint(2,false));
    data.push_back(DataPoint(3,true));
    std::ostringstream dataStream;
    writer.writeData(dataStream, data, ",");
    ASSERT_EQUAL(dataStream.str(), "1,3\n", "StatisticsWriter::writeData failed to skip inactive data points");

    DataVector disabledData;
    disabledData.push_back(DataPoint(1,true));
    StatisticsTrackerOptions disabledOptions;
    disabledOptions.doCov = disabledOptions.doMax = disabledOptions.doMean = disabledOptions.doMin = disabledOptions.doVar = true;
    MultivariateTracker disabledTracker(1, disabledOptions);
    disabledTracker.update(disabledData);

    std::ostringstream disabledHistograms;
    writer.writeHistograms(disabledHistograms, disabledTracker);
    ASSERT(disabledHistograms.str().find("No valid data") != std::string::npos, "StatisticsWriter::writeHistograms failed to report no data for a disabled histogram");

}

/**
 * Test StatisticsApp file input and output
 */
void test_statisticsapp_fileio()
{

    std::string fileInput = "test_statisticsapp_fileio_input.csv";
    std::string fileOutput = "test_statisticsapp_fileio_output.csv";

    std::ofstream input(fileInput.c_str());
    input << "1,2" << std::endl;
    input << "3,4" << std::endl;
    input.close();

    CommandLineOptions options;
    options.fileInput = fileInput;
    options.fileOutput = fileOutput;
    options.delimiter = ",";
    options.showStatistics = true;

    StatisticsApp app(options);
    int status = app.run();
    ASSERT_EQUAL(status, 0, "StatisticsApp::run failed on file input");
    app.display();

    std::ifstream output(fileOutput.c_str());
    std::stringstream contents;
    contents << output.rdbuf();
    output.close();

    ASSERT(contents.str().find("Statistics") != std::string::npos, "StatisticsApp::display failed to write statistics to output file");

    std::remove(fileInput.c_str());
    std::remove(fileOutput.c_str());

}

/**
 * Test StatisticsApp skips blank lines and comment rows
 */
void test_statisticsapp_blankcomment()
{

    std::string fileInput = "test_statisticsapp_blankcomment_input.csv";
    std::string fileOutput = "test_statisticsapp_blankcomment_output.csv";

    std::ofstream input(fileInput.c_str());
    input << "# comment row" << std::endl;
    input << "1,2" << std::endl;
    input << std::endl;
    input << "3,4" << std::endl;
    input.close();

    CommandLineOptions options;
    options.fileInput = fileInput;
    options.fileOutput = fileOutput;
    options.delimiter = ",";
    options.comment = "#";
    options.showStatistics = true;

    StatisticsApp app(options);
    int status = app.run();
    ASSERT_EQUAL(status, 0, "StatisticsApp::run failed on blank line and comment row input");
    app.display();

    std::ifstream output(fileOutput.c_str());
    std::stringstream contents;
    contents << output.rdbuf();
    output.close();

    ASSERT(contents.str().find("2") != std::string::npos, "StatisticsApp::run failed to process rows around blank line and comment");

    std::remove(fileInput.c_str());
    std::remove(fileOutput.c_str());

}

/**
 * Test Logger::log at ERROR, WARNING, and INFO levels
 */
void test_logger_levels()
{

    Logger::Level::Type saved = Logger::logLevel;
    Logger::logLevel = Logger::Level::INFO;

    std::ostringstream error;
    std::streambuf * oldErrBuf = std::cerr.rdbuf(error.rdbuf());
    Logger::log(Logger::Level::ERROR) << "error message" << std::endl;
    std::cerr.rdbuf(oldErrBuf);
    ASSERT(error.str().find("ERROR") != std::string::npos, "Logger::log failed to label ERROR level");

    std::ostringstream warning;
    std::streambuf * oldWarnBuf = std::cerr.rdbuf(warning.rdbuf());
    Logger::log(Logger::Level::WARNING) << "warning message" << std::endl;
    std::cerr.rdbuf(oldWarnBuf);
    ASSERT(warning.str().find("WARNING") != std::string::npos, "Logger::log failed to label WARNING level");

    std::ostringstream info;
    std::streambuf * oldInfoBuf = std::cout.rdbuf(info.rdbuf());
    Logger::log(Logger::Level::INFO) << "info message" << std::endl;
    std::cout.rdbuf(oldInfoBuf);
    ASSERT(info.str().find("INFO") != std::string::npos, "Logger::log failed to label INFO level");

    Logger::logLevel = saved;

}


void test_statisticsapp_badoutputfile()
{

    std::string fileInput = "test_statisticsapp_badoutputfile_input.csv";
    std::ofstream input(fileInput.c_str());
    input << "1,2" << std::endl;
    input.close();

    CommandLineOptions options;
    options.fileInput = fileInput;
    options.fileOutput = "/no_such_directory_xyz/out.csv";
    options.delimiter = ",";
    options.showStatistics = true;

    std::ostringstream fallback;
    std::streambuf * oldCoutBuf = std::cout.rdbuf(fallback.rdbuf());
    StatisticsApp app(options);
    app.run();
    app.display();
    std::cout.rdbuf(oldCoutBuf);

    ASSERT(fallback.str().find("Statistics") != std::string::npos, "StatisticsApp failed to fall back to standard output on bad output file");

    std::remove(fileInput.c_str());

}

/**
 * Test StatisticsApp exits immediately on blank line when blankEOF is set
 */
void test_statisticsapp_blankeof()
{

    std::string fileInput = "test_statisticsapp_blankeof_input.csv";
    std::string fileOutput = "test_statisticsapp_blankeof_output.csv";

    std::ofstream input(fileInput.c_str());
    input << "1,2" << std::endl;
    input << std::endl;
    input << "3,4" << std::endl;
    input.close();

    CommandLineOptions options;
    options.fileInput = fileInput;
    options.fileOutput = fileOutput;
    options.delimiter = ",";
    options.blankEOF = true;
    options.showStatistics = true;

    StatisticsApp app(options);
    app.run();
    app.display();

    std::ifstream output(fileOutput.c_str());
    std::stringstream contents;
    contents << output.rdbuf();
    output.close();

    ASSERT(contents.str().find("Count") != std::string::npos, "StatisticsApp::run failed to produce statistics before blank EOF line");
    ASSERT(contents.str().find("3") == std::string::npos, "StatisticsApp::run failed to stop processing at blank EOF line");

    std::remove(fileInput.c_str());
    std::remove(fileOutput.c_str());

}

/**
 * Test StatisticsApp skip/keep line counts and uniform row sampling
 */
void test_statisticsapp_rowselection()
{

    std::string fileInputSkipKeep = "test_statisticsapp_rowselection_skipkeep_input.csv";
    std::string fileOutputSkipKeep = "test_statisticsapp_rowselection_skipkeep_output.csv";

    std::ofstream inputSkipKeep(fileInputSkipKeep.c_str());
    inputSkipKeep << "1,1" << std::endl;
    inputSkipKeep << "2,2" << std::endl;
    inputSkipKeep << "3,3" << std::endl;
    inputSkipKeep << "4,4" << std::endl;
    inputSkipKeep << "5,5" << std::endl;
    inputSkipKeep.close();

    CommandLineOptions optionsSkipKeep;
    optionsSkipKeep.fileInput = fileInputSkipKeep;
    optionsSkipKeep.fileOutput = fileOutputSkipKeep;
    optionsSkipKeep.delimiter = ",";
    optionsSkipKeep.numLinesToSkip = 1;
    optionsSkipKeep.numLinesToKeep = 2;
    optionsSkipKeep.showStatistics = true;
    optionsSkipKeep.statisticsOptions.doCov = optionsSkipKeep.statisticsOptions.doMax = optionsSkipKeep.statisticsOptions.doMean = optionsSkipKeep.statisticsOptions.doMin = optionsSkipKeep.statisticsOptions.doVar = true;

    StatisticsApp appSkipKeep(optionsSkipKeep);
    appSkipKeep.run();
    appSkipKeep.display();

    std::ifstream outputSkipKeep(fileOutputSkipKeep.c_str());
    std::stringstream contentsSkipKeep;
    contentsSkipKeep << outputSkipKeep.rdbuf();
    outputSkipKeep.close();

    std::string outputContentSkipKeep = contentsSkipKeep.str();
    ASSERT(outputContentSkipKeep.find("2.000000") != std::string::npos, "StatisticsApp::run failed to skip leading lines");
    ASSERT(outputContentSkipKeep.find("5.000000") == std::string::npos, "StatisticsApp::run failed to stop after keep count reached");

    std::remove(fileInputSkipKeep.c_str());
    std::remove(fileOutputSkipKeep.c_str());

    std::string fileInputSampling = "test_statisticsapp_rowselection_sampling_input.csv";
    std::string fileOutputSampling = "test_statisticsapp_rowselection_sampling_output.csv";

    std::ofstream inputSampling(fileInputSampling.c_str());
    inputSampling << "1,1" << std::endl;
    inputSampling << "2,2" << std::endl;
    inputSampling << "3,3" << std::endl;
    inputSampling << "4,4" << std::endl;
    inputSampling.close();

    CommandLineOptions optionsSampling;
    optionsSampling.fileInput = fileInputSampling;
    optionsSampling.fileOutput = fileOutputSampling;
    optionsSampling.delimiter = ",";
    optionsSampling.sampling.mode = SamplerOptions::UNIFORM;
    optionsSampling.sampling.step = 2;
    optionsSampling.showStatistics = true;

    StatisticsApp appSampling(optionsSampling);
    appSampling.run();
    appSampling.display();

    std::ifstream outputSampling(fileOutputSampling.c_str());
    std::stringstream contentsSampling;
    contentsSampling << outputSampling.rdbuf();
    outputSampling.close();

    ASSERT(contentsSampling.str().find("2") != std::string::npos, "StatisticsApp::run failed to process rows with uniform sampling");

    std::remove(fileInputSampling.c_str());
    std::remove(fileOutputSampling.c_str());

}

/**
 * Test StatisticsApp applies header row titles and falls back to numbered titles on mismatch
 */
void test_statisticsapp_headertitles()
{

    std::string fileInputMatch = "test_statisticsapp_headertitles_match_input.csv";
    std::string fileOutputMatch = "test_statisticsapp_headertitles_match_output.csv";

    std::ofstream inputMatch(fileInputMatch.c_str());
    inputMatch << "col1,col2" << std::endl;
    inputMatch << "1,2" << std::endl;
    inputMatch << "3,4" << std::endl;
    inputMatch.close();

    CommandLineOptions optionsMatch;
    optionsMatch.fileInput = fileInputMatch;
    optionsMatch.fileOutput = fileOutputMatch;
    optionsMatch.delimiter = ",";
    optionsMatch.headerRow = 1;
    optionsMatch.showStatistics = true;

    StatisticsApp appMatch(optionsMatch);
    appMatch.run();
    appMatch.display();

    std::ifstream outputMatch(fileOutputMatch.c_str());
    std::stringstream contentsMatch;
    contentsMatch << outputMatch.rdbuf();
    outputMatch.close();

    ASSERT(contentsMatch.str().find("col1") != std::string::npos, "StatisticsApp::run failed to apply header row titles");

    std::remove(fileInputMatch.c_str());
    std::remove(fileOutputMatch.c_str());

    std::string fileInputMismatch = "test_statisticsapp_headertitles_mismatch_input.csv";
    std::string fileOutputMismatch = "test_statisticsapp_headertitles_mismatch_output.csv";

    std::ofstream inputMismatch(fileInputMismatch.c_str());
    inputMismatch << "col1,col2,col3" << std::endl;
    inputMismatch << "1,2" << std::endl;
    inputMismatch << "3,4" << std::endl;
    inputMismatch.close();

    CommandLineOptions optionsMismatch;
    optionsMismatch.fileInput = fileInputMismatch;
    optionsMismatch.fileOutput = fileOutputMismatch;
    optionsMismatch.delimiter = ",";
    optionsMismatch.headerRow = 1;
    optionsMismatch.showStatistics = true;

    StatisticsApp appMismatch(optionsMismatch);
    appMismatch.run();
    appMismatch.display();

    std::ifstream outputMismatch(fileOutputMismatch.c_str());
    std::stringstream contentsMismatch;
    contentsMismatch << outputMismatch.rdbuf();
    outputMismatch.close();

    ASSERT(contentsMismatch.str().find("col1") == std::string::npos, "StatisticsApp::run failed to discard mismatched header titles");

    std::remove(fileInputMismatch.c_str());
    std::remove(fileOutputMismatch.c_str());

}

/**
 * Test StatisticsApp strict-parsing dimension mismatch skipping and string/numeric filter rejection
 */
void test_statisticsapp_rowfiltering()
{

    std::string fileInputStrict = "test_statisticsapp_rowfiltering_strict_input.csv";
    std::string fileOutputStrict = "test_statisticsapp_rowfiltering_strict_output.csv";

    std::ofstream inputStrict(fileInputStrict.c_str());
    inputStrict << "1,2" << std::endl;
    inputStrict << "3,4,5" << std::endl;
    inputStrict << "6,7" << std::endl;
    inputStrict.close();

    CommandLineOptions optionsStrict;
    optionsStrict.fileInput = fileInputStrict;
    optionsStrict.fileOutput = fileOutputStrict;
    optionsStrict.delimiter = ",";
    optionsStrict.strictParsing = true;
    optionsStrict.showStatistics = true;
    optionsStrict.statisticsOptions.doCov = optionsStrict.statisticsOptions.doMax = optionsStrict.statisticsOptions.doMean = optionsStrict.statisticsOptions.doMin = optionsStrict.statisticsOptions.doVar = true;

    StatisticsApp appStrict(optionsStrict);
    appStrict.run();
    appStrict.display();

    std::ifstream outputStrict(fileOutputStrict.c_str());
    std::stringstream contentsStrict;
    contentsStrict << outputStrict.rdbuf();
    outputStrict.close();

    ASSERT(contentsStrict.str().find("2.000000") != std::string::npos, "StatisticsApp::run failed to process rows around a strict dimension mismatch");

    std::remove(fileInputStrict.c_str());
    std::remove(fileOutputStrict.c_str());

    std::string fileInputFilter = "test_statisticsapp_rowfiltering_filter_input.csv";
    std::string fileOutputFilter = "test_statisticsapp_rowfiltering_filter_output.csv";

    std::ofstream inputFilter(fileInputFilter.c_str());
    inputFilter << "foo,1" << std::endl;
    inputFilter << "bar,100" << std::endl;
    inputFilter << "foo,100" << std::endl;
    inputFilter << "foo,2" << std::endl;
    inputFilter.close();

    CommandLineOptions optionsFilter;
    optionsFilter.fileInput = fileInputFilter;
    optionsFilter.fileOutput = fileOutputFilter;
    optionsFilter.delimiter = ",";
    optionsFilter.filterRows.addStringFilter(0,"bar",true,true,false);
    optionsFilter.filterRows.addNumericFilter(1,0,10,true);
    optionsFilter.showStatistics = true;
    optionsFilter.statisticsOptions.doCov = optionsFilter.statisticsOptions.doMax = optionsFilter.statisticsOptions.doMean = optionsFilter.statisticsOptions.doMin = optionsFilter.statisticsOptions.doVar = true;

    StatisticsApp appFilter(optionsFilter);
    appFilter.run();
    appFilter.display();

    std::ifstream outputFilter(fileOutputFilter.c_str());
    std::stringstream contentsFilter;
    contentsFilter << outputFilter.rdbuf();
    outputFilter.close();

    ASSERT(contentsFilter.str().find("1.500000") != std::string::npos, "StatisticsApp::run failed to process rows passing string and numeric filters");

    std::remove(fileInputFilter.c_str());
    std::remove(fileOutputFilter.c_str());

}

/**
 * Test StatisticsApp writes raw filtered data and reads carriage-return line feed terminated lines
 */
void test_statisticsapp_io()
{

    std::string fileInputRaw = "test_statisticsapp_io_raw_input.csv";
    std::string fileOutputRaw = "test_statisticsapp_io_raw_output.csv";

    std::ofstream inputRaw(fileInputRaw.c_str());
    inputRaw << "1,2" << std::endl;
    inputRaw.close();

    CommandLineOptions optionsRaw;
    optionsRaw.fileInput = fileInputRaw;
    optionsRaw.fileOutput = fileOutputRaw;
    optionsRaw.delimiter = ",";
    optionsRaw.showFilteredData = true;

    StatisticsApp appRaw(optionsRaw);
    appRaw.run();
    appRaw.display();

    std::ifstream outputRaw(fileOutputRaw.c_str());
    std::stringstream contentsRaw;
    contentsRaw << outputRaw.rdbuf();
    outputRaw.close();

    ASSERT_EQUAL(contentsRaw.str(), "1 2\n", "StatisticsApp::run failed to write raw filtered data");

    std::remove(fileInputRaw.c_str());
    std::remove(fileOutputRaw.c_str());

    std::string fileInputCrlf = "test_statisticsapp_io_crlf_input.csv";
    std::string fileOutputCrlf = "test_statisticsapp_io_crlf_output.csv";

    std::ofstream inputCrlf(fileInputCrlf.c_str(), std::ios::binary);
    inputCrlf << "1,2\r\n";
    inputCrlf << "3,4\r\n";
    inputCrlf.close();

    CommandLineOptions optionsCrlf;
    optionsCrlf.fileInput = fileInputCrlf;
    optionsCrlf.fileOutput = fileOutputCrlf;
    optionsCrlf.delimiter = ",";
    optionsCrlf.showStatistics = true;
    optionsCrlf.statisticsOptions.doCov = optionsCrlf.statisticsOptions.doMax = optionsCrlf.statisticsOptions.doMean = optionsCrlf.statisticsOptions.doMin = optionsCrlf.statisticsOptions.doVar = true;

    StatisticsApp appCrlf(optionsCrlf);
    appCrlf.run();
    appCrlf.display();

    std::ifstream outputCrlf(fileOutputCrlf.c_str());
    std::stringstream contentsCrlf;
    contentsCrlf << outputCrlf.rdbuf();
    outputCrlf.close();

    ASSERT(contentsCrlf.str().find("2.000000") != std::string::npos, "StatisticsApp::run failed to read carriage-return line feed terminated lines");

    std::remove(fileInputCrlf.c_str());
    std::remove(fileOutputCrlf.c_str());

}

/**
 * Test StatisticsApp displays covariance, correlation, least squares, and histogram outputs
 */
void test_statisticsapp_alloutputs()
{

    std::string fileInput = "test_statisticsapp_alloutputs_input.csv";
    std::string fileOutput = "test_statisticsapp_alloutputs_output.csv";

    std::ofstream input(fileInput.c_str());
    input << "1,2" << std::endl;
    input << "2,4" << std::endl;
    input << "3,6" << std::endl;
    input.close();

    CommandLineOptions options;
    options.fileInput = fileInput;
    options.fileOutput = fileOutput;
    options.delimiter = ",";
    options.showCovariance = true;
    options.showCorrelation = true;
    options.showLeastSquaresOffset = true;
    options.showLeastSquaresSlope = true;
    options.showHistogram = true;
    options.statisticsOptions.doCov = options.statisticsOptions.doMax = options.statisticsOptions.doMean = options.statisticsOptions.doMin = options.statisticsOptions.doVar = true;
    options.statisticsOptions.histogramOptions.enabled = true;
    options.statisticsOptions.histogramOptions.binCount = 2;
    options.statisticsOptions.histogramOptions.cacheSize = 2;

    StatisticsApp app(options);
    app.run();
    app.display();

    std::ifstream output(fileOutput.c_str());
    std::stringstream contents;
    contents << output.rdbuf();
    output.close();

    std::string outputContent = contents.str();
    ASSERT(outputContent.find("Covariance") != std::string::npos, "StatisticsApp::display failed to write covariance matrix");
    ASSERT(outputContent.find("Correlation") != std::string::npos, "StatisticsApp::display failed to write correlation matrix");
    ASSERT(outputContent.find("Offset") != std::string::npos, "StatisticsApp::display failed to write least squares offset matrix");
    ASSERT(outputContent.find("Slope") != std::string::npos, "StatisticsApp::display failed to write least squares slope matrix");
    ASSERT(outputContent.find("Histogram") != std::string::npos, "StatisticsApp::display failed to write histograms");

    std::remove(fileInput.c_str());
    std::remove(fileOutput.c_str());

}

/**
 * Unit tester for clistats
 */
int
main()
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
    tests.push_back(test_dynamichistogram_degenerate);
    tests.push_back(test_dynamichistogram_expand_after_degenerate);
    tests.push_back(test_datavector);
    tests.push_back(test_datafilters);
    tests.push_back(test_statisticstracker);
    tests.push_back(test_multivariatetracker);
    tests.push_back(test_multivariatetracker_inactive);
    tests.push_back(test_multivariatetracker_zero_variance_correlation);
    tests.push_back(test_stringsplitter_tointegers_filter_edges);
    tests.push_back(test_rowsampler_uniform);
    tests.push_back(test_rowsampler_random);
    tests.push_back(test_statisticsapp_fileio);
    tests.push_back(test_statisticsapp_blankcomment);
    tests.push_back(test_commandlineparser_flags_io);
    tests.push_back(test_commandlineparser_flags_filters);
    tests.push_back(test_commandlineparser_flags_sampling);
    tests.push_back(test_commandlineparser_flags_output);
    tests.push_back(test_commandlineparser_help_version);
    tests.push_back(test_commandlineparser_invalid);
    tests.push_back(test_pdrv);
    tests.push_back(test_dynamichistogram_order);
    tests.push_back(test_statisticswriter_output);
    tests.push_back(test_logger_levels);
    tests.push_back(test_statisticsapp_badoutputfile);
    tests.push_back(test_statisticsapp_blankeof);
    tests.push_back(test_statisticsapp_rowselection);
    tests.push_back(test_statisticsapp_headertitles);
    tests.push_back(test_statisticsapp_rowfiltering);
    tests.push_back(test_statisticsapp_io);
    tests.push_back(test_statisticsapp_alloutputs);

    // ======================================================================
    // Define test metrics
    // ----------------------------------------------------------------------
    unsigned long numTestsTotal = static_cast<unsigned long>(tests.size());
    unsigned long numTestsPass = 0;
    unsigned long numTestsFail = 0;

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
    std::cout << "        Testing " << ApplicationProperties::NAME  << " " << ApplicationProperties::version() << std::endl;
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
        // LCOV_EXCL_START
        catch (std::exception const & ex)
        {
            numTestsFail++;
            errors.push_back(ex.what());
            message = "F";
        }
        // LCOV_EXCL_STOP
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
    // LCOV_EXCL_START
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
    // LCOV_EXCL_STOP
}
