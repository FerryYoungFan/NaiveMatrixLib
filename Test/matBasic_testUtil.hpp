#pragma once

#include <iostream>
#include <chrono>


class TestTimer
{
public:
    TestTimer() = default;
    ~TestTimer() = default;

    void tic()
    {
        m_running = true;
        m_time_start = std::chrono::steady_clock::now();
    }

    double toc(const char *printInfo = nullptr)
    {
        double res = 0.0;
        if (m_running)
        {
            m_time_end = std::chrono::steady_clock::now();
            m_running = false;
            auto time_diff = m_time_end - m_time_start;
            auto microSecs = std::chrono::duration_cast<std::chrono::microseconds>(time_diff).count();
            res = static_cast<double>(microSecs) / 1.0e6;
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
    std::chrono::steady_clock::time_point m_time_start;
    std::chrono::steady_clock::time_point m_time_end;

    bool m_running{false};
};
