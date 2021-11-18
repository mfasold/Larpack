/**
 * @file DataCollector.hpp Contains a simple class to collect information during the applied methods
 * @author $Author: Jan $ 
 * @author Jan Bruecker
 * @date $Date: 2008-06-13 14:58:41 +0200 (Fri, 13 Jun 2008) $
 *
 */
#ifndef _DATACOLLECTOR_
#define _DATACOLLECTOR_

#include <string>
#include <map>
#include <boost/any.hpp>
#include <boost/noncopyable.hpp>
#include <boost/pool/detail/singleton.hpp>

namespace larrpack {
  
  typedef std::map<std::string,  boost::any> DataOutputTable;

  /**
   * @class DataCollector
   * @brief Global (Singleton) class that collects data that is needed 
   *        globally and may be written to an output file
   * 
   */
  class DataCollector : private boost::noncopyable // Make sure we have only a single DataCollector instance
  {
  public:
    static DataCollector& instance();
      
    DataCollector();
    ~DataCollector();

    /// Returns data entry if value is of type T (throws boost::bad_any_cast exception otherwise)
    template<class T>
    T& getAs(const std::string key) {
      return boost::any_cast<T&>(mDataTable[key]);
    }

    // @note The following function does not allow for string casting,
    //       hence we need multiple overloaded insert functions
    //     void insert(std::string key, const boost::any& value)
    //     {
    //       mDataTable[key] = value;
    //     }

      
    
    void insert(std::string key, double value); // Insert an IntensityType
    void insert(std::string key, std::string value); // Insert a string
    void insert(std::string key, size_t value); // Insert a size_t
    void insert(std::string key, int value); // Insert a int

    /// Insert an existing table of string/IntensityType pairs.
    void insert(std::map<std::string, double> statistics);


    /// Write the collected Data to a file
    void writeToFile(const char*  filename);

    /// Empty complete logger
    void clear();
    

    //       std::string toString();
    //       void insert(std::string, IntensityType value);
    //       void insert(DataOutputTable addTable);
    //       void write(const char *  filename);
  private:
    DataOutputTable mDataTable;
  };

}
#endif
