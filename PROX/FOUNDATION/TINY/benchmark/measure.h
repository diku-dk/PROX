#ifndef MEASURE_H
#define MEASURE_H

#include <util_log.h>

#include <map>
#include <vector>
#include <string>

class Measure
  {
      using Samples = std::map<size_t, double>;
      using Rows = std::map<std::string, Samples>;
      using Tables = std::map<std::string, Rows>;

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
      auto begin = measures.begin();
      auto end = measures.end();

      for (auto i = begin; i != end; ++i)
      {
        typedef Rows::iterator rowiterator;
        auto rbegin = (*i).second.begin();
        auto rend = (*i).second.end();

        for (auto j = rbegin; j != rend; ++j)
        {
          logging << (*i).first << " " << (*j).first << "\t";

          for (size_t k = 1; k <= (*j).second.size() ; ++k)
          {
            typedef Samples::iterator siterator;
            auto s = (*j).second.find(k);
            logging << (*s).second << "\t";
          }
          logging << util::Log::newline();
        }
      }
    }
  };

// MEASURE_H
#endif
