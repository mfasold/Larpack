/**
 * @file DataCollector.cpp Defines methods to collect different kind of data in a global (singleton) hash 
 * @author $Author: Jan $ 
 * @author Jan Bruecker
 * @date 
 *
 */
#include <iostream>
#include <fstream>

#include <map>
#include <string>
#include <iterator> 
#include "HookcurveStatistics.hpp"
#include "DataCollector.hpp"
#include <boost/any.hpp>
#include <boost/pool/detail/singleton.hpp>

using namespace std;
using namespace larrpack;
using namespace boost;

/**
 * Constructor initiallizes mDataTable, so it's available before the start of main.
 * 
 */
DataCollector::DataCollector() : mDataTable()
{
}

DataCollector::~DataCollector()
{}

/**
 * Accessor to the only instance of the class
 * 
 */
DataCollector& DataCollector::instance() 
{
  return boost::details::pool::singleton_default<DataCollector>::instance();
}

// /**
//  * Returns the entry of the map given the key (only if the value is a string, otherwise an exception is thrown)
//  * 
//  * @param key
//  * @return value of the given key
//  * 
//  */
// string DataCollector::getString(string key)
// {
//   string back = any_cast<string> (mDataTable[key]);
//   return back;
// }


/**
 * Creates a new map entry for a pair of a key (string) and value (IntensityType)
 * 
 * @param key Id or key of the value to be inserted
 * @param value Value of data entry to be inserted
 */
void DataCollector::insert(string key, IntensityType value)
{
  mDataTable[key] = (boost::any)value;
}

/**
 * Creates a new map entry for a pair of a key (string) and value (string)
 * 
 * @param key Id or key of the value to be inserted
 * @param value Value of data entry to be inserted
 */
void DataCollector::insert(string key, string value)
{
  mDataTable[key] =   (boost::any)value;
}

/**
 * Creates a new map entry for a pair of a key (string) and value (size_t)
 * 
 * @param key Id or key of the value to be inserted
 * @param value Value of data entry to be inserted
 */
void DataCollector::insert(string key, size_t value)
{
  mDataTable[key] =   (boost::any)value;
}

/**
 * Creates a new map entry for a pair of a key (string) and value (int)
 * 
 * @param key Id or key of the value to be inserted
 * @param value Value of data entry to be inserted
 */
void DataCollector::insert(string key, int value)
{
  mDataTable[key] =   (boost::any)value;
}


/**
 * Inserts the map of string/IntensityType pairs into the data table
 * 
 * @param statistics
 * 
 */
void DataCollector::insert(IntensityTypeStatistics statistics)
{
  mDataTable.insert(statistics.begin(), statistics.end());
}

/**
 * Write the data table to file in the format  "key\tvalue\n"
 * 
 * @param filename
 * 
 */
void DataCollector::writeToFile(const char* filename)
{
  ofstream logfile;
  logfile.open(filename);
  for (DataOutputTable::const_iterator itr = mDataTable.begin(); 
       itr != mDataTable.end(); ++itr) {
    // Write key
    logfile << itr->first << "\t";

    // Write value depending on variable type
    const boost::any&  value = mDataTable[itr->first];
    const std::type_info& type= value.type();
 
    if (type == typeid(std::string))
      logfile << any_cast<string>(value);
    else if (type == typeid(double))
      logfile << any_cast<double>(value);
    else if (type == typeid(int))
      logfile << any_cast<int>(value);
    else if (type == typeid(size_t))
      logfile << any_cast<size_t>(value);
    else
      logfile << "Value type is not considered yet";
    logfile << endl;
  }
  logfile.close();
}

/**
 * Removes all elements
 * 
 */
void DataCollector::clear() {
  mDataTable.clear();  
}
