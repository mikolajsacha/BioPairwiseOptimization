#include <mutex>
#include <condition_variable>
#include "../include/semaphore.h"

Semaphore::Semaphore(int count) : count(count) {}
Semaphore::Semaphore() : count(0) {}

void Semaphore::notify() {
    std::unique_lock<std::mutex> lck(mtx);
    ++count;
    cv.notify_one();
}

void Semaphore::wait() {
    std::unique_lock<std::mutex> lck(mtx);
    while (count == 0) {
          cv.wait(lck);
    }
    --count;
}
