#pragma once

#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>

class PersistentWorker {
private:
    std::thread worker_thread;
    std::mutex mtx;
    std::condition_variable cv_main;
    std::condition_variable cv_worker;

    bool terminate = false;
    bool work_ready = false;
    bool work_done = true;

    std::function<void()> current_task;

public:
    PersistentWorker() {
        worker_thread = std::thread([this]() {
            while (true) {
                std::unique_lock<std::mutex> lock(mtx);
                // Sleep until there is work or the program is shutting down
                cv_worker.wait(lock, [this]() { return work_ready || terminate; });

                if (terminate) break;

                // Execute
                if (current_task) {
                    current_task();
                }

                // Mark as finished and wake up the main thread
                work_ready = false;
                work_done = true;
                lock.unlock();
                cv_main.notify_one();
            }
        });
    }

    ~PersistentWorker() {
        {
            std::lock_guard<std::mutex> lock(mtx);
            terminate = true;
        }
        cv_worker.notify_one();
        if (worker_thread.joinable()) {
            worker_thread.join();
        }
    }

    void launch(std::function<void()> task) {
        std::lock_guard<std::mutex> lock(mtx);
        current_task = task;
        work_ready = true;
        work_done = false;
        cv_worker.notify_one(); // Wake up the background thread
    }

    void wait() {
        std::unique_lock<std::mutex> lock(mtx);
        // Block the main thread until the worker signals completion
        cv_main.wait(lock, [this]() { return work_done; });
    }
};
