#include <mutex>
#include <condition_variable>

/*!< Custom implementation of a simple Semaphore using mutex and contition_variable */
class Semaphore
{
    private:
        std::mutex mtx;
        std::condition_variable cv;
        int count;

    public:
        Semaphore(int count);
        Semaphore();

        void notify();
        void wait();
};
