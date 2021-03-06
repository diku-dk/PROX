#ifndef UTIL_TIMER_H
#define UTIL_TIMER_H

#include <chrono>

#include <cassert>


namespace util 
{

  class Timer
  {
  protected:

    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_end;

    bool m_start_called;

  public:

    Timer()
    : m_start_called(false)
    {}

  public:

    void start()
    {
      m_end = m_start = std::chrono::high_resolution_clock::now();
      m_start_called = true;
    }

    void stop()
    {
      assert( m_start_called || !"stop(): start was never invokved");

      m_end = std::chrono::high_resolution_clock::now();

      assert( m_end >= m_start || !"stop(): std::chrono::high_resolution_clock::now() must have returned a bad-value?");

      m_start_called = false;
    }

    float duration()
    {
      assert( m_start_called || !"stop(): start was never invokved");

      this->stop();

      assert( m_end >= m_start || !"duration(): end time was smaller than start time");

      return std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start).count();
    }

    /**
     * Get Time.
     *
     * @return     The time between start and stop in milliseconds.
     */
    float operator()()const
    {
      assert( m_end >= m_start || !"operator(): Must call stop before operator()");

      return std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start).count();
    }
  };

}  // end of namespace util

// UTIL_TIMER_H
#endif 
