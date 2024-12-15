//
//  logger_V2.cpp
//  Final_Project_V3
//
//  Created by jonathon hope on 12/14/24.
//

#include "logger_V2.hpp"

// Initialize static members
Logger *Logger::instancePtr = nullptr;
mutex Logger::mtx;

Logger *Logger::getInstance(const string &filename)
{
    if (nullptr == instancePtr)
    {
        lock_guard<mutex> lock(mtx);
        if (nullptr == instancePtr)
        {
            instancePtr = new Logger(filename.empty() ? "logger.log" : filename);
        }
    }
    return instancePtr;
}

Logger::Logger(const string &filename)
{
    openLogFile(filename);
}

Logger::~Logger()
{
    logFile.close();
}

void Logger::setLogFile(const string &filename)
{
    if (nullptr == instancePtr)
    {
        getInstance(filename);
        return;
    }
    logFile.close();
    openLogFile(filename);
}

void Logger::openLogFile(const string &filename)
{
    logFile.open(filename, ios::app);
    if (!logFile.is_open())
    {
        cerr << "error opening log file" << endl;
    }
}

void Logger::log(LogLevel level, const string &message)
{
    //const string message_output = message.str();
    lock_guard<mutex> lock(mtx);

    // Get current timestamp
    time_t now = time(0);
    tm *timeinfo = localtime(&now);
    char timestamp[20];
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", timeinfo);

    ostringstream logEntry;
    if (level == FUNCTION_START || level == FUNCTION_END || level == PROGRAM_START || level == PROGRAM_END)
    {
        // Create log entry
        logEntry << "[" << timestamp << "] " << levelToString(level) << ": " << message << endl;
        cout << "[" << timestamp << "] " << levelToString(level) << ": " << message << endl;

    }
    
    else if(level == BRANCH_COUNT){
        cout << "[" << timestamp << "] " << levelToString(level) << ": " << message << endl;
    }
    else
    {
        logEntry << levelToString(level) << message << endl;
        cout << levelToString(level) << message << endl;

    }

    // Output to console (can be disabled with `#define USE_LOGGING_OUT`)
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

string Logger::levelToString(LogLevel level)
{
    switch (level)
    {
    case HIT_POINCARE:
        return "Number of branches that reached the poincare section at  ";
    case HIT_MOON:
        return "Number of branches that didn't hit the Moon: ";
    case POS_VS_NEG:
        return ""; // You may want to output specific data here
    case MEASUREMENTS:
        return ""; // You may want to output specific data here
    case SC_STATE:
        return "";
    case PO_DATA:
        return "";
    case MAX_DELTA_V:
        return "";
    case TRANSFER_INFO:
        return "";
        case BRANCH_COUNT:
            return "Currently on branch number: ";
    case FUNCTION_START:
        return "Function: ";
    case FUNCTION_END:
        return "Function: ";
    case PROGRAM_START:
        return "*********Program Starting*********";
    case PROGRAM_END:
        return "*********Program END*********";
    default:
        return "UNKNOWN";
    }
}
