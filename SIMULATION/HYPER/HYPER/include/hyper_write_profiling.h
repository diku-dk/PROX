#ifndef HYPER_WRITE_PROFILING_H
#define HYPER_WRITE_PROFILING_H

#include <util_profiling.h>
#include <util_matlab_write_profiling.h>
#include <util_log.h>

#include <string>
#include <vector>
#include <fstream>   // needed for std::ofstream

namespace hyper
{

  bool write_profiling(std::string const & filename)
  {
    util::Log        logging;

    std::string const newline = util::Log::newline();

    std::ofstream matlab;

    matlab.open(filename.c_str(),std::ios::out);

    if(! matlab.is_open())
    {
      logging << "hyper::write_profiling(): error could not open file = " << filename.c_str() << util::Log::newline();

      return false;
    }

    matlab << "close all;" << std::endl;
    matlab << "clear all;" << std::endl;
    matlab << "clc;"       << std::endl;

    matlab << util::matlab_write_profiling()  << std::endl;

    // 20XX-YY-ZZ Robert: These times are kernel times (if OpenCL is used).
    // They are part of narrow_time. Subtract this sum from narrow_time to
    // obtain the amount of time spent on pre-/post-processing and data transfer.

    matlab << "contact_time =  tandem_traversal - exact_test;";

    //--- setting up 'eyecandy' part -----------------------------------------
    matlab <<  "red   = [0.7,0.1,0.1];"   << std::endl;
    matlab <<  "green = [0.1,0.7,0.1];" << std::endl;
    matlab <<  "blue  = [0.1,0.1,0.7];"  << std::endl;
    matlab <<  "f_type = 'Times';"      << std::endl;
    matlab <<  "f_size = 20;"           << std::endl;

    //--- plotting -----------------------------------------------------------

    matlab << "figure(300);" << std::endl;
    matlab << "plot(contacts);" << std::endl;
    matlab << "title('Contacts', 'FontSize', f_size, 'FontName', f_type); " << std::endl;
    matlab << "xlabel('Frame', 'FontSize', f_size, 'FontName', f_type); " << std::endl;
    matlab << "ylabel('Number of contacts', 'FontSize', f_size, 'FontName', f_type)" << std::endl;
    matlab << "print(gcf,'-depsc2','contacts');" << std::endl;
    matlab << "print(gcf,'-dpng','contacts');" << std::endl;

    matlab.flush();
    matlab.close();

    logging << "hyper::write_profiling(): Done writing profile data..." << newline;

    return true;
  }


} // namespace hyper

// HYPER_WRITE_PROFILING_H
#endif
