/*
 * File:   bb_sqlite.h
 * Author: karol
 *
 * Created on 4 listopad 2009, 17:43
 */

#ifndef _BB_SQLITE_H
#define _BB_SQLITE_H

#include "std.h"
#include "serialsaveload.h"
#include "exceptions.h"
#include "boost.h"
#include "sqlite3.h"
#include "serialhash.h"
#include "proxies.h"

namespace boostbase
{

/**
 * @brief SQLite database for binary serializable objects
 *
 * With this class one can transparently store and retrieve objects
 * (instances of arbitrary classes) using an SQLite database stored
 * on disk. Any class can be stored, given that it supplies
 * the operator|(serializer &, object &), as defined in serializer.h
 *
 * The objects are stored in individual files, compressed with zlib
 * in a separate directory, identified by date,time and an MD5 hash.
 * The SQLite database stores metadata as needed.
 *
 */


class base
{
    ///SQLite database handle
    sqlite3 * sqlite_db;
    ///log
    std::stringstream string_log;
    std::ostream * sqlite_log;
    ///directory for storing files
    fs::path base_dir;

    ///thread handles for writing to files
    std::list<std::thread>  thread_pool;
    ///mutexes
    std::mutex dbaccess;
    std::mutex logaccess;

    /// evaluates a non-select command 'statement' such as INSERT
    void command(const std::string & statement)
    {

        std::lock_guard<std::mutex> lock(dbaccess);

        log() << "command(): " << statement << std::endl;
        int result = 0;
        sqlite3_stmt * stmt;

        while(true)
        {
            result = sqlite3_prepare_v2(sqlite_db, statement.c_str(), statement.size(), &stmt, 0);
            if(result == SQLITE_BUSY || result == SQLITE_LOCKED)
                SLEEP_MS(20);
            else break;
        }

        if(result != SQLITE_OK && result != SQLITE_DONE)
        {
            log() << sqlite3_errmsg(sqlite_db) << std::endl;
            throw exception::sqlite_command_error();
        }

        while(true)
        {
            //#//pragma omp critical
            result = sqlite3_step(stmt);
            if(result == SQLITE_BUSY || result == SQLITE_LOCKED)
                SLEEP_MS(20);
            else break;
        }

        if(result != SQLITE_OK && result != SQLITE_DONE)
        {
            log() << sqlite3_errmsg(sqlite_db) << std::endl;
            throw exception::sqlite_command_error();
        }

        while(true)
        {
            sqlite3_finalize(stmt);
            if(result == SQLITE_BUSY || result == SQLITE_LOCKED)
                SLEEP_MS(20);
            else break;
        }
    }

    /// evaluates a SELECT statement 'statement', given 'columns' number of columns
    /// returns a vector of rows in the form of vectors of strings
    std::vector<std::vector<std::string> > eval_select(const std::string & statement, int columns)
    {
        std::lock_guard<std::mutex> lock(dbaccess);

        log() << "eval_select(): " << statement << std::endl;
        std::vector<std::vector<std::string> > items;
        int result = 0;
        sqlite3_stmt * stmt;

        while(true)
        {
            result = sqlite3_prepare_v2(sqlite_db, statement.c_str(), statement.size(), &stmt, 0);
            if(result == SQLITE_BUSY || result == SQLITE_LOCKED)
                SLEEP_MS(20);
            else break;
        }

        if(result != SQLITE_OK && result != SQLITE_DONE)
        {
            log() << sqlite3_errmsg(sqlite_db) << std::endl;
            throw exception::sqlite_command_error();
        }

        while(result != SQLITE_DONE)
        {
            std::vector<std::string> cur_row;

            while(true)
            {
                //#//pragma omp critical
                result = sqlite3_step(stmt);
                if(result == SQLITE_BUSY || result == SQLITE_LOCKED)
                    SLEEP_MS(20);
                else break;
            }

            if(result == SQLITE_ROW)
            {
                for(int i = 0; i < columns; i++)
                {
                    cur_row.push_back((const char*) sqlite3_column_text(stmt, i));
                    log() << cur_row.back() << ", ";
                }
                log() << std::endl;
            }
            if(cur_row.size() != columns)
            {
                //std::cerr << "eval_select(): should have selected " << columns << " columns, got " << cur_row.size() << std::endl;
            }
            else
                items.push_back(cur_row);
        }
        while(true)
        {
            sqlite3_finalize(stmt);
            if(result == SQLITE_BUSY || result == SQLITE_LOCKED)
                SLEEP_MS(20);
            else break;
        }
        return items;
    }

    /// returns a vector of objects matching 'query'
    /// checks the object's MD5 sum before including it in the results
    /// NOTE: you need to be certain of the stored objects' type, otherwise
    /// you will get garbage data
    template<class item_t>
    std::vector<item_t> eval_get(const std::string & query)
    {
        std::vector<std::vector<std::string> > selected;
        try
        {
            selected = eval_select(query, 2);
        }
        catch(exception::sqlite_command_error & e)
        {
            log() << "get(): select couldn't be evaluated\n";
            return std::vector<item_t>();
        }
        std::vector<item_t> results;

        for(std::vector<std::string> & field : selected)
        {
            if(field.size() == 2)
                try
                {
                    results.push_back(serialload<item_t > (fs::path(base_dir / field[0]), field[1]));
                }
                catch(exception::file_not_found e)
                {
                    log() << "get(): object file not found: " << field[0] << std::endl;
                }
                catch(exception::wrong_md5 e)
                {
                    log() << "get(): md5 mismatch: " << field[1] << std::endl;
                }
        }
        log() << "get(): restored " << results.size() << " objects\n";
        return results;
    }

public:

    /// constructor
    base(const fs::path & dbfile, const fs::path & dir, bool readonly = false) :
        base_dir(dir)
    {
        sqlite_log = &string_log;

        //critical section
        {
            std::lock_guard<std::mutex> lock(dbaccess);

            int result = 0;
            if(!dbfile.parent_path().empty() && !fs::exists(dbfile.parent_path()))
                fs::create_directories(dbfile.parent_path());
            if(!fs::exists(base_dir))
                fs::create_directories(base_dir);


            //do {
            if(readonly)
                result = sqlite3_open_v2(dbfile.string().c_str(), &sqlite_db, SQLITE_OPEN_READONLY, "unix-none");
            else
                result = sqlite3_open_v2(dbfile.string().c_str(), &sqlite_db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, "unix-none");
            //} while(result==SQLITE_BUSY || result==SQLITE_LOCKED);

            if(result != SQLITE_OK && result != SQLITE_DONE)
            {
                log() << sqlite3_errmsg(sqlite_db) << std::endl;
                sqlite3_close(sqlite_db);
                throw exception::sqlite_failed_opening_db();
            }
        }
        //end of critical section

        try
        {
            command("CREATE TABLE IF NOT EXISTS objects (id_o INTEGER PRIMARY KEY AUTOINCREMENT, filename TEXT, date DATETIME, md5 TEXT)");
        }
        catch(exception::sqlite_command_error e)
        {
            log() << e.what() << std::endl;
            throw exception::sqlite_failed_creating_db();
        }

    }

    /// the destructor cares for clearing the thread pool
    /// and ensures all objects are stored before exiting
    ~base()
    {
        for(std::thread & t : thread_pool)
            if(t.joinable())
                t.join();
        thread_pool.clear();
        sqlite3_close(sqlite_db);
    }

    /// retrieve objects without writing SQL directly
    /// usage:
    /// get<item_t>(where(a,b)(c,d)(e,f));
    /// will return a vector of items for which
    /// a=b,c=d and e=f
    template<class item_t>
    std::vector<item_t> get(const std::vector<pair_t> & wheres)
    {
        std::stringstream s;
        s << "SELECT filename,md5 from objects\nWHERE ";
        for(const pair_t & p : wheres)
        {
            std::vector<std::string> typeandvalue;
            boost::split(typeandvalue, p.second, boost::is_any_of(":"));

            s << p.first << "=\'" << typeandvalue[1] << "\'";
            if(p != wheres.back())
                s << " AND ";
        }
        std::vector<item_t> results;
        //#//pragma omp critical
        results = eval_get<item_t>(s.str());
        return results;
    }


    /// retrieve objects without writing SQL directly
    /// usage:
    /// get<item_t>(between(a,b,c)(d,e,f))
    /// will return a vector of items for which
    /// a >= b and a < c and d >= e and d < f
    ///
    template<class item_t>
    std::vector<item_t> get(const std::vector<tween_t> & betweens)
    {
        std::stringstream s;
        s << "SELECT filename,md5 from objects\nWHERE ";
        for(const tween_t & p : betweens)
        {
            s << p.field << " BETWEEN \'" << p.begin << "\' AND \'" << p.end << "\'";
            if(p != betweens.back())
                s << " AND ";
        }
        std::vector<item_t> results;
        //#//pragma omp critical
        results = eval_get<item_t>(s.str());
        return results;
    }

    /// retrieve objects without writing SQL directly
    /// see above functions
    ///
    template<class item_t>
    std::vector<item_t> get(const std::vector<pair_t> & wheres, const std::vector<tween_t> & betweens)
    {
        std::stringstream s;
        s << "SELECT filename,md5 from objects\nWHERE ";
        for(const pair_t & p : wheres)
        {
            std::vector<std::string> typeandvalue;
            boost::split(typeandvalue, p.second, boost::is_any_of(":"));

            s << p.first << "=\'" << typeandvalue[1] << "\'";
            if(p != wheres.back())
                s << " AND ";
        }
        if(betweens.size() != 0) s << " AND ";
        for(const tween_t & p : betweens)
        {
            s << p.field << " BETWEEN \'" << p.begin << "\' AND \'" << p.end << "\'";
            if(p != betweens.back())
                s << " AND ";
        }
        std::vector<item_t> results;
        //#//pragma omp critical
        results = eval_get<item_t>(s.str());
        return results;
    }

    /// retrieves only metadata ('columns') for which the
    /// WHERE clause specified by 'wheres' matches,
    /// in the form of a vector of textual rows
    std::vector<std::vector<std::string> > get_metadata(const std::vector<std::string> & columns, const std::vector<pair_t> & wheres)
    {
        std::stringstream s;
        s << "SELECT ";
        for(int i = 0; i < columns.size(); i++)
        {
            const std::string & c = columns[i];
            s << c ;
            if(i != columns.size() - 1)
                s << ",";
        }
        s << " from objects\n";

        if(wheres.size())
            s << "WHERE ";

        for(const pair_t & p : wheres)
        {
            std::vector<std::string> typeandvalue;
            boost::split(typeandvalue, p.second, boost::is_any_of(":"));

            s << p.first << "=\'" << typeandvalue[1] << "\'";
            if(p != wheres.back())
                s << " AND ";
        }

        return eval_select(s.str(), columns.size());
    }

    /// retrieves only metadata ('columns') for which both the
    /// WHERE and BETWEEN clauses specified by 'wheres' and 'betweens' match,
    /// in the form of a vector of textual rows
    std::vector<std::vector<std::string> > get_metadata(const std::vector<std::string> & columns, const std::vector<pair_t> & wheres, const std::vector<tween_t> & betweens)
    {
        std::stringstream s;
        s << "SELECT ";
        for(int i = 0; i < columns.size(); i++)
        {
            const std::string & c = columns[i];
            s << c ;
            if(i != columns.size() - 1)
                s << ",";
        }
        s << " from objects\n";

        if(wheres.size())
            s << "WHERE ";
        for(const pair_t & p : wheres)
        {
            std::vector<std::string> typeandvalue;
            boost::split(typeandvalue, p.second, boost::is_any_of(":"));

            s << p.first << "=\'" << typeandvalue[1] << "\'";
            if(p != wheres.back())
                s << " AND ";
        }
        if(betweens.size())
            s << " AND ";
        for(const tween_t & p : betweens)
        {
            s << p.field << " BETWEEN \'" << p.begin << "\' AND \'" << p.end << "\'";
            if(p != betweens.back())
                s << " AND ";
        }
        return eval_select(s.str(), columns.size());
    }

    /// retrieves only metadata ('columns') for which the
    /// BETWEEN clause specified by 'betweens' matches,
    /// in the form of a vector of textual rows
    std::vector<std::vector<std::string> > get_metadata(const std::vector<std::string> & columns,  const std::vector<tween_t> & betweens)
    {
        std::stringstream s;
        s << "SELECT ";
        for(int i = 0; i < columns.size(); i++)
        {
            const std::string & c = columns[i];
            s << c ;
            if(i != columns.size() - 1)
                s << ",";
        }
        s << " from objects\n";

        if(betweens.size())
            s << "WHERE ";

        for(const tween_t & p : betweens)
        {
            s << p.field << " BETWEEN \'" << p.begin << "\' AND \'" << p.end << "\'";
            if(p != betweens.back())
                s << " AND ";
        }

        return eval_select(s.str(), columns.size());
    }

    /// retrieves a vector of column names in the object table
    std::vector<std::string> columns()
    {
        std::vector<std::vector<std::string> > qry = eval_select(std::string("pragma table_info(objects);"), 3);
        std::vector<std::string> res;
        for(std::vector<std::string> & v : qry)
        {
            res.push_back(v[1]);
        }
        return res;
    }

    /// deletes objects from the database (not from disk) for which
    /// the query specified by 'wheres' and 'betweens' matches
    void remove(const std::vector<pair_t> & wheres, const std::vector<tween_t> & betweens)
    {
        std::stringstream s;
        s << "DELETE from objects\nWHERE ";
        for(const pair_t  & p : wheres)
        {
            std::vector<std::string> typeandvalue;
            boost::split(typeandvalue, p.second, boost::is_any_of(":"));

            s << p.first << "=\'" << typeandvalue[1] << "\'";
            if(p != wheres.back())
                s << " AND ";
        }
        if(betweens.size() != 0) s << " AND ";
        for(const tween_t & p :  betweens)
        {
            s << p.field << " BETWEEN \'" << p.begin << "\' AND \'" << p.end << "\'";
            if(p != betweens.back())
                s << " AND ";
        }
        command(s.str());
    }

    /// the blocking version of the save function, avoid
    template<class item_t>
    void store_thread(const item_t & item, const std::vector<pair_t> & column_values)
    {

        log() << "store_thread id = " << std::this_thread::get_id() << ", this = " << reinterpret_cast<void *>(this) << std::endl;
        std::stringstream s1, s2; // dwie części polecenia SQL, w jednej nazwy kolumn, w drugiej wartości
        std::stringstream f; // nazwa pliku
        std::string md5 = md5gen<item_t> (item); // unikalny hash obiektu
        pt::ptime local_time = pt::microsec_clock::local_time(); // aktualny czas z dokładnością do us
        std::string fn = pt::to_simple_string(local_time) + "_" + md5 + ".z";

        //f << base_dir << "/";

        // tworzymy zapytanie
        s1 << "INSERT INTO objects ( ";
        s2 << " VALUES ( ";

        // sortujemy kolumny alfabetycznie
        std::vector<pair_t> sorted_column_values = column_values;
        std::sort(boost::begin(sorted_column_values), boost::end(sorted_column_values), [](const pair_t &a, const pair_t &b) { return a.first <= b.first; });

        std::vector<std::string> existing_columns = columns();

        for(pair_t & p : sorted_column_values)
        {
            // uzupełniamy niesitniejące kolumny
            std::vector<std::string>    typeandvalue;
            boost::split(typeandvalue, p.second, boost::is_any_of(":"));

            if(std::find(existing_columns.begin(),existing_columns.end(),p.first) == existing_columns.end())
            {
                //nie znaleziono
                std::stringstream create;
                create << "ALTER TABLE objects ADD COLUMN " << p.first << " ";
                if(typeandvalue[0] == typeid(std::string).name())
                    create << "TEXT";
                if(typeandvalue[0] == typeid(double).name())
                    create << "REAL";
                if(typeandvalue[0] == typeid(int).name())
                    create << "REAL";
                if(typeandvalue[0] == typeid(float).name())
                    create << "REAL";

                try
                {
                    command(create.str());
                }
                catch(exception::sqlite_command_error & e)
                {
                    log() << e.what() << std::endl;
                }
            }

            f << p.first << "=" << typeandvalue[1] << ",";
            s1 << p.first << ',';
            s2 << '\'' << typeandvalue[1] << "\' ,";
        }
        s1 << "filename, date, md5 ) ";
        f << "md5=" << fn;

        s2 << '\'' << /*f.str()*/ fn << "\' , \'" << pt::to_simple_string(local_time) << "\', \'" << md5 << "\' )";

        // wykonanie polecenia SQL
        // #// pragma omp critical
        try
        {
            command((s1.str() + s2.str()));
        }
        catch(exception::sqlite_command_error & e)
        {
            log() << e.what() << std::endl;
        }

        // zapisanie obiektu do pliku
        serialsave<item_t> (item, fs::path(base_dir / /*f.str()*/fn));
    }

    /// stores object 'item' in the database, with the metadata by which it can
    /// be later identified specified in 'column_values', which is a vector of pairs
    /// of std::string (understood as key = value)
    ///
    /// this function executes the inserion by spawning a thread and exiting immediately,
    /// deferring the actual I/O to a separate thread, allowing the application to
    /// continue without blocking
    template<class item_t>
    void store(const item_t & item, const std::vector<pair_t> & column_values)
    {
        /// if curious about the syntax, it's a lambda function, i.e. a function
        /// object defined in place, which is then passed as the function to execute
        /// in the new thread, which is then pushed to the thread pool, so we can
        /// retain its handle and allow it to safely terminate
        thread_pool.push_back(std::thread::thread(
                                  [=]{ store_thread<item_t>(item,column_values); }
        ));
    }

    /// in/out function for accessing the SQLite log
    std::ostream & log()
    {
        std::lock_guard<std::mutex> lock(logaccess);
        return *sqlite_log;
    }
    /// specify a different stream for logging (e.g. &std::cout)
    void SetStream(std::ostream * o)
    {
        sqlite_log = o;
    }
};

}



#endif  /* _BB_SQLITE_H */

