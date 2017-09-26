#include <util/threadpool.h>
//#include <util/string.h>

#include <iostream>









//------------------------------------------------------------------------------
// Task
//------------------------------------------------------------------------------

/*!
	Destructor.
 */
Task::~Task()
{
}

//------------------------------------------------------------------------------
// ThreadPool::Worker
//------------------------------------------------------------------------------

/*!
	Thread pool worker class.

	Workers represent threads, each processing it's assigned tasks.
 */
class ThreadPool::Worker {
public:
	Worker(ThreadPool *pool);
	void run();

private:
	ThreadPool *m_pool;     ///< Thread pool the worker belongs to
};

/*!
	Constructor.

	\param pool Thread pool this worker belongs to
 */
ThreadPool::Worker::Worker(ThreadPool *pool) :
	m_pool(pool)
{
}

/*!
	Main loop for the worker threads.

	Waits for a task, runs it, repeats.
 */
void ThreadPool::Worker::run()
{
	Task *task;

	while (true) {
		{
			std::unique_lock<std::mutex> lock(m_pool->m_queue_mutex);

			// Wait for task
			while (!m_pool->m_stop && m_pool->m_tasks.empty())
				m_pool->m_start_cond.wait(lock);

			// Exit if thread pool is stopped
			if (m_pool->m_stop)
				return;

			// Get the task from the queue
			task = m_pool->m_tasks.front();
			m_pool->m_tasks.pop_front();
			m_pool->m_activeTasks++;
		}

		// Run task
		task->run();

		// Check if finished
		{
			std::lock_guard<std::mutex> lock(m_pool->m_queue_mutex);

			m_pool->m_activeTasks--;
			if (m_pool->m_activeTasks == 0)
				m_pool->m_finish_cond.notify_one();
		}
	}
}

//------------------------------------------------------------------------------
// ThreadPool
//------------------------------------------------------------------------------

static void startWorker(void *arg)
{
	static_cast<ThreadPool::Worker *>(arg)->run();
}

/*!
	Constructor.

	\param numWorkers Number of workers to use (0 = use number of system cores)
 */
ThreadPool::ThreadPool(int numWorkers) :
	m_activeTasks(0),
	m_stop(false)
{
	if (numWorkers <= 0)
	{
		numWorkers = getNumSystemCores();
	}

	//char* pPath;
	//pPath = getenv ("OMP_NUM_THREADS");
	//if( pPath )
	//	numWorkers = fromString<int>(std::string(pPath));

	if (singleThreaded)
		numWorkers = 1;

	std::cout << "Threadpool: working with " << numWorkers << " threads\n";

	for (int i = 0; i < numWorkers; i++)
		m_workers.push_back(Worker(this));

	for (int i = 0; i < numWorkers; i++)
		m_threads.push_back(new std::thread(startWorker, &m_workers[i]));
}

/*!
	Destructor.
 */
ThreadPool::~ThreadPool()
{
	m_stop = true;
	m_start_cond.notify_all();

	for (unsigned i = 0; i < m_threads.size(); i++)
	{
		m_threads[i]->join();
		delete m_threads[i];
	}
}

/*!
	Enqueues a single task to be processed.

	\param task Task to process
 */
void ThreadPool::enqueue(Task *task)
{
	// add task
	std::lock_guard<std::mutex> lock(m_queue_mutex);
	m_tasks.push_back(task);

	// wake up one worker
	m_start_cond.notify_one();
}

/*!
	Waits (blocks) until all enqueued tasks have been processed.
 */
void ThreadPool::waitForAll()
{
	std::unique_lock<std::mutex> lock(m_queue_mutex);

	while (m_activeTasks > 0 || !m_tasks.empty())
		m_finish_cond.wait(lock);
}

/*!
	\return Returns the number of workers.
 */
const int ThreadPool::getNumWorkers() const
{
	return m_workers.size();
}

/*!
	\return Returns the number of available system cores.
 */
int ThreadPool::getNumSystemCores()
{
	return std::thread::hardware_concurrency();
}


bool ThreadPool::singleThreaded = false;


ThreadPool *g_instance = NULL;  ///< Global thread pool

/*!
	\return Returns the global thread pool.
 */
ThreadPool* ThreadPool::instance()
{
	if (!g_instance)
		g_instance = new ThreadPool();
	return g_instance;
}



