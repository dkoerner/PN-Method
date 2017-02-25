#pragma once

#include <iostream>
#include <vector>
#include <deque>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <sstream>
#include <util/timer.h>

#include <math/RNG.h>

class ThreadPool;

//! Abstract task base class.
/*!
	Tasks are run by worker threads and should implement the run() method to
	implement their work.
 */
class Task {
public:
	virtual ~Task();

	/*!
		Called to perform the task's work.
	 */
	virtual void run() = 0;
};

//! Thread pool.
/*!
	Implements a simple thread pool where a set of tasks is processed by
	a number of individual worker threads.
 */
class ThreadPool {
public:
	class Worker;

	ThreadPool(int numWorkers = 0);
	~ThreadPool();

	void enqueue(Task *task);

	template <class T>
	void enqueue(const std::vector<T *> &tasks);

	template <class T>
	void enqueueAndWait(const std::vector<T *> &tasks);

	void waitForAll();

	const int getNumWorkers() const;

	static int getNumSystemCores();

	static ThreadPool* instance();

	static bool singleThreaded;

private:
	friend class ThreadPool::Worker;

	std::vector<Worker> m_workers;                  ///< List of workers
	std::vector<std::thread *> m_threads;       ///< List of threads

	std::deque<Task *> m_tasks;                     ///< Queue of tasks
	int m_activeTasks;                              ///< Number of active tasks

	std::mutex m_queue_mutex;                   ///< Mutex to protect the queue
	std::condition_variable m_start_cond;       ///< Start condition variable
	std::condition_variable m_finish_cond;      ///< Stop condition variable
	bool m_stop;
};

/*!
	Enqueues a set of tasks to be processed.

	\param tasks Tasks to process
 */
template <class T>
void ThreadPool::enqueue(const std::vector<T *> &tasks)
{
	// Add task
	std::lock_guard<std::mutex> lock(m_queue_mutex);

	m_tasks.insert(m_tasks.end(), tasks.begin(), tasks.end());

	// Wake up all workers
	m_start_cond.notify_all();
}

/*!
	Enqueues a set of tasks and waits (blocks) until they have all been
	processed.

	\param tasks Tasks to process and wait for
 */
template <class T>
void ThreadPool::enqueueAndWait(const std::vector<T *> &tasks)
{
	enqueue(tasks);
	waitForAll();
}




// allows to switch between number of samples or given computation time (for equal time comparisons)
struct Terminator
{
	enum EMode
	{
		ETime,
		ENumSamples
	};

	Terminator( int numSamples ):
		m_mode(ENumSamples),
		m_numSamples(numSamples),
		m_sample(0)
	{
	}

	Terminator( double seconds ):
		m_mode(ETime),
		m_seconds(seconds)
	{
	}

	bool keepRunning()
	{
		if( m_mode == ENumSamples )
			return m_sample < m_numSamples;
		return m_timer.elapsedSeconds() < m_seconds;
	}

	void advance()
	{
		++m_sample;
	}

	void reset()
	{
		m_sample = 0;
		m_timer.reset();
		m_timer.start();
	}

	void done()
	{
		m_timer.stop();
	}

	void printReport()
	{
		std::cout << "Took " << m_timer.elapsedSeconds() << "s for " << m_sample << " samples." << std::endl;
	}


	EMode m_mode;

	// run till number of samples
	int m_numSamples;
	int m_sample;

	// run till time limit is reached
	Timer m_timer;
	double m_seconds;
};


struct MonteCarloTaskInfo
{
	MonteCarloTaskInfo(int taskid, int numTasks):
		taskid(taskid),
		numTasks(numTasks),
		samples(0),
		rng(123+taskid)
	{
	}

	static MonteCarloTaskInfo create( int taskid, int numTasks )
	{
		return MonteCarloTaskInfo(taskid, numTasks);
	}

	// per task info ---
	int taskid;
	int numTasks;
	int samples;
	RNGd rng;
};


template<typename TaskInfo>
struct GenericTask : Task
{
	typedef std::function<void(TaskInfo&)> RunMethod;

	GenericTask( TaskInfo& info, RunMethod _run )
		: Task(),
		  m_info(info),
		  m_run(_run)
	{
	}
	virtual ~GenericTask()
	{
	}

	virtual void run()
	{
		m_run( m_info );
	}


	TaskInfo m_info;
	RunMethod m_run;
};

template<typename TaskInfo>
void runGenericTasks( typename GenericTask<TaskInfo>::RunMethod run, Terminator& term, int numThreads)
{
	std::cout << "runTasks" << std::endl;
	//int numThreads = ThreadPool::getNumSystemCores();

	std::vector<GenericTask<TaskInfo>*> tasks;
	for (int i = 0; i < numThreads; i++)
		tasks.push_back(new GenericTask<TaskInfo>(TaskInfo::create(i, numThreads), run));

	for( term.reset();term.keepRunning();term.advance() )
	{
		//std::cout << term.m_sample <<" " << term.m_numSamples <<  std::endl;
		if(numThreads>1)
			ThreadPool::instance()->enqueueAndWait(tasks);
		else
			tasks[0]->run();
	}
	term.done();
	term.printReport();

	for( auto task:tasks )
		delete task;
}

/*
template<typename TaskInfo>
void runGenericTasks( typename GenericTask<TaskInfo>::RunMethod run, int numThreads)
{
	std::cout << "runTasks" << std::endl;
	//int numThreads = ThreadPool::getNumSystemCores();

	std::vector<GenericTask<TaskInfo>*> tasks;
	for (int i = 0; i < numThreads; i++)
		tasks.push_back(new GenericTask<TaskInfo>(TaskInfo::create(i, numThreads), run));

	if(numThreads>1)
		ThreadPool::instance()->enqueueAndWait(tasks);
	else
		tasks[0]->run();

	for( auto task:tasks )
		delete task;
}
*/
