#pragma once

#include <iostream>
#include <sys/time.h> // Performance test timer


class TestTimer
{
public:
    TestTimer() = default;
    ~TestTimer() = default;

    void tic()
    {
        m_running = true;
        gettimeofday(&m_time_start, NULL);
    }

    double toc(const char *printInfo = nullptr)
    {
        double res = 0.0;
        if (m_running)
        {
            gettimeofday(&m_time_end, NULL);
            m_running = false;
            res = static_cast<double>(m_time_end.tv_sec - m_time_start.tv_sec) + static_cast<double>(m_time_end.tv_usec - m_time_start.tv_usec) / 1000000.0;
            if (nullptr != printInfo)
            {
                std::cout << "Time elapsed - " << printInfo << ": ";
                if (res > 1.0)
                {
                    std::cout << res << " s\n";
                }
                else if (res > 0.001)
                {
                    std::cout << res * 1000.0 << " ms\n";
                }
                else
                {
                    std::cout << res * 1000000.0 << " us\n";
                }
            }
        }
        return res;
    }

private:
    struct timeval m_time_start
    {
        0, 0
    };

    struct timeval m_time_end
    {
        0, 0
    };

    bool m_running{false};
};
