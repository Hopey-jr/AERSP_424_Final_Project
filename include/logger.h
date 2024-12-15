//
//  logger.h
//  Final_Project_V3
//
//  Created by jonathon hope on 12/14/24.
//

#pragma once
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <mutex>
using namespace std;

// Enum to represent log levels
enum LogLevel
{
    HIT_POINCARE,
    HIT_MOON,
    STABLE_VS_UNSTABLE,
    MEASUREMENTS,
    SC_STATE,
    PO_DATA,
    MAX_DELTA_V,
    TRANSFER_INFO,
    FUNCTION_START,
    FUNCTION_END
};

class Logger
{
private:
    // Static pointer to the Singleton instance
    static Logger *instancePtr;

    // Mutex to ensure thread safety
    static mutex mtx;

public:
    // Deleting the copy constructor to prevent copies
    Logger(const Logger &obj) = delete;

    // Static method to get the Singleton instance
  
    
    static Logger *getInstance(const string &filename = "")
    {
        if (nullptr == instancePtr)
        {
            lock_guard<mutex> lock(mtx);
            if (nullptr == instancePtr)
            {
                instancePtr = new Logger(filename.empty() ? "logger.log" : filename.c_str());
            }
        }
        return instancePtr;
    }

    // Constructor: Opens the log file in append mode
    Logger(const string &filename)
    {
        openLogFile(filename);
    }

    // Destructor: Closes the log file
    ~Logger() { logFile.close(); }

    void setLogFile(const string &filename)
    {
        if (nullptr == instancePtr)
        {
            getInstance(filename);
            return;
        }
        logFile.close();
        openLogFile(filename);
    }

    // Logs a message with a given log level
    void log(LogLevel level, const ostringstream  &message)
    {
        const string message_output = message.str();
        lock_guard<mutex> lock(mtx);
        // Get current timestamp
        time_t now = time(0);
        tm *timeinfo = localtime(&now);
        char timestamp[20];
        strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", timeinfo);

        ostringstream logEntry;
        if(level == FUNCTION_START || level == FUNCTION_END){
            // Create log entry
            
            logEntry << "[" << timestamp << "] " << levelToString(level) << ": " << message_output << endl;
        }
        else{
           logEntry << levelToString(level) << message_output << endl;
        }
        // Output to console
#ifdef USE_LOGGING_OUT
        cout << logEntry.str();
#endif

        // Output to log file
        if (logFile.is_open())
        {
            logFile << logEntry.str();
            logFile.flush(); // Ensure immediate write to file
        }
    }

private:
    void openLogFile(const string &filename)
    {
        logFile.open(filename, ios::app);
        if (!logFile.is_open())
        {
            cerr << "error opening log file" << endl;
        }
    }

    ofstream logFile; // File stream for the log file

    // Converts log level to a string for output
    string levelToString(LogLevel level)
    {
        switch (level)
        {
        case HIT_POINCARE:
            return "Number of branches that hit the manifold: ";
        case HIT_MOON:
            return "Number of branches that hit the Moon: ";
        case STABLE_VS_UNSTABLE:
            return ""; //I want the output to be within the sentence
        case MEASUREMENTS:
            return ""; //I want to output multiple pieces of data
        case SC_STATE:
            return "";
        case PO_DATA:
            return "";
        case MAX_DELTA_V:
                return "";
        case TRANSFER_INFO:
                return "";
        case FUNCTION_START:
                return "Function: ";
        case FUNCTION_END:
                return "Function: ";
        default:
            return "UNKNOWN";
        }
    }
};

// Initialize static members
Logger *Logger::instancePtr = nullptr;
mutex Logger::mtx;
