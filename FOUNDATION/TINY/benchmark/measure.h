#ifndef MEASURE_H
#define MEASURE_H

#include <util_log.h>

#include <map>
#include <vector>
#include <string>

class Measure
  {
    typedef std::map<size_t,double>       Samples;
    typedef std::map<std::string,Samples> Rows;
    typedef std::map<std::string,Rows>    Tables;
    
    std::string table;
    std::string row;
    
  private:
    
    Tables measures;
  
  public:
    
    void SelectOperation(std::string op)
    {
      row = op;
    }
    
    void SelectType(std::string type)
    {
      table = type;
    }
    
    double& operator [] (const size_t idx)
    {
      return measures[table][row][idx];
    }
    
    void clear()
    {
      measures.clear();
    }
    
    void dump()
    {
      util::Log logging;

      typedef Tables::iterator iterator;
      iterator begin = measures.begin();
      iterator end = measures.end();
      
      for (iterator i = begin ; i!=end; ++i)
      {
        typedef Rows::iterator rowiterator;
        rowiterator rbegin = (*i).second.begin();
        rowiterator rend = (*i).second.end();
        
        for (rowiterator j = rbegin ; j!=rend; ++j)
        {  
          logging << (*i).first << " " << (*j).first << "\t";
          
          for (size_t k = 1; k <= (*j).second.size() ; ++k)
          {
            typedef Samples::iterator siterator;
            siterator s = (*j).second.find(k); 
            logging << (*s).second << "\t";
          }
          logging << util::Log::newline();
        }
      }
    }
  };

// MEASURE_H
#endif
