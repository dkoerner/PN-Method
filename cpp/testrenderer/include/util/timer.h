#pragma once


#include <chrono>
#include <ostream>
#include <type_traits>

////std::this_thread::sleep_for(std::chrono::seconds(2));

struct Timer
{
	typedef std::conditional< std::chrono::high_resolution_clock::is_steady,
	std::chrono::high_resolution_clock,
	std::chrono::steady_clock >::type clock_type ;
	typedef std::chrono::milliseconds milliseconds;
	typedef std::chrono::duration<float> seconds;

	explicit Timer(bool run = false):
		m_running(false),
		m_elapsed_seconds(0.0f)
	{
		if(run)
			start();
	}
	void reset()
	{
		m_running = false;
		m_elapsed_seconds = 0.0f;
	}
	void start()
	{
		m_start_time = clock_type::now();
		m_running = true;
	}
	void stop()
	{
		m_elapsed_seconds += std::chrono::duration_cast<seconds>(clock_type::now() - m_start_time).count();
		m_running = false;
	}

//	milliseconds elapsedMiliseconds() const
//	{
//		return std::chrono::duration_cast<milliseconds>(clock_type::now() - start_time);
//	}

	float elapsedSeconds()
	{
		if(m_running)
		{
			stop();
			start();
		}
		return m_elapsed_seconds;
	}
//	template <typename T, typename Traits>
//	friend std::basic_ostream<T, Traits>& operator<<(std::basic_ostream<T, Traits>& out, const Timer& timer)
//	{
//		return out << timer.elapsedSeconds();
//	}
private:
	clock_type::time_point m_start_time;
	float                  m_elapsed_seconds;
	bool                   m_running;
};
