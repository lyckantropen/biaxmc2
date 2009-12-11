/* 
 * File:   bb_sqlite.h
 * Author: karol
 *
 * Created on 4 listopad 2009, 17:43
 */

#ifndef _BB_SQLITE_H
#define	_BB_SQLITE_H

#include "std.h"
#include "serialsaveload.h"
#include "exceptions.h"
#include "boost.h"
#include "sqlite3.h"
#include "serialhash.h"
#include "proxies.h"

namespace boostbase {

    /**
     * @brief Interfejs SQL
     *
     * Klasa ta służy do przechowywania obiektów na dysku i odtwarzania ich z bazy danych.
     * By to osiągnąć, serializuje obiekty, pakuje je za pomocą zlib i zapisuje informacje
     * w bazie danych SQLite. Obiekty zapisywane są na dysku w osobnych plikach. Obiekty
     * są rozróżniane na podstawie wzorca MD5.
     *
     * Aby skorzystać z klasy, trzeba zaimplementować zewnętrznie
     * operator|(serializer &, obiekt &)
     */

    class base {
        ///obiekt bazy SQL
        sqlite3 * sqlite_db;
        ///wiadomości wewnętrzne SQLite, można się do nich dostać przez log()
        std::stringstream sqlite_log;
        ///bazowy katalog bazy danych
        fs::path base_dir;

        /**
         * @brief ewaluacja komendy innej niż zapytanie
         * @param statement komenda SQL
         */
        void command(const std::string & statement) {
            sqlite_log << "command(): " << statement << std::endl;
            int result = 0;
            sqlite3_stmt * stmt;

	    int k=0;
            do {
		k++;
                result = sqlite3_prepare_v2(sqlite_db, statement.c_str(), statement.size(), &stmt, 0);
		if(k%10000==0){
                    std::cout << "Waiting for database lock 1\n";
                    std::cout << sqlite3_errmsg(sqlite_db) << std::endl;
                }
            } while(result==SQLITE_BUSY || result==SQLITE_LOCKED);

            if (result != SQLITE_OK) {
                sqlite_log << sqlite3_errmsg(sqlite_db) << std::endl;
                std::cout << sqlite_log.str() << std::endl;
                throw exception::sqlite_command_error();
            }

	    k=0;
            do {
		k++;
                result = sqlite3_step(stmt);
		if(k%1000==0)
			std::cout << "Waiting for database lock 2\n";
            } while(result==SQLITE_BUSY || result==SQLITE_LOCKED);

            if (result == SQLITE_ERROR) {
                sqlite_log << sqlite3_errmsg(sqlite_db) << std::endl;
                std::cout << sqlite_log.str() << std::endl;
                throw exception::sqlite_command_error();
            }

	k=0;
            do {
		k++;
                sqlite3_finalize(stmt);
		if(k%1000==0)
			std::cout << "Waiting for database lock 3\n";
            } while(result==SQLITE_BUSY || result==SQLITE_LOCKED);
        }

        /**
         * @brief ewaluacja zapytania SQL
         * @param statement zapytanie
         * @param columns liczba kolumn w zapytaniu
         * @return tabela wyników w postaci tekstowej
         */
        std::vector<std::vector<std::string> > eval_select(const std::string & statement, int columns) {
            sqlite_log << "eval_select(): " << statement << std::endl;
            std::vector<std::vector<std::string> > items;
            int result = 0;
            sqlite3_stmt * stmt;

            do {
                result = sqlite3_prepare_v2(sqlite_db, statement.c_str(), statement.size(), &stmt, 0);
            } while(result==SQLITE_BUSY || result==SQLITE_LOCKED);

            if (result != SQLITE_OK) {
                sqlite_log << sqlite3_errmsg(sqlite_db) << std::endl;
                //std::cout << sqlite_log.str() << std::endl;
                throw exception::sqlite_command_error();
            }

            while (result != SQLITE_DONE) {
                std::vector<std::string> cur_row;

                do {
                    result = sqlite3_step(stmt);
                } while(result==SQLITE_BUSY || result==SQLITE_LOCKED);

                if (result == SQLITE_ROW) {
                    for (int i = 0; i < columns; i++) {
                        cur_row.push_back((const char*) sqlite3_column_text(stmt, i));
                        sqlite_log << cur_row.back() << ", ";
                    }
                    sqlite_log << std::endl;
                }
                if (cur_row.size() != columns) {
                    sqlite_log << "eval_select(): should have selected " << columns << " columns, got " << cur_row.size() << std::endl;
                } else
                    items.push_back(cur_row);
            }
            do {
                sqlite3_finalize(stmt);
            } while(result==SQLITE_BUSY || result==SQLITE_LOCKED);
            
            return items;
        }

        /**
         * @brief ewaluacja zapytania, zwracająca konkretne obiekty
         * @param query zapytanie SQL
         * @return wektor obiektów zwróconych przez zapytanie
         */
        template<class item_t>
        std::vector<item_t> eval_get(const std::string & query){
            std::vector<std::vector<std::string> > selected;
            selected = eval_select(query, 2);
            std::vector<item_t> results;

            foreach(std::vector<std::string> & field, selected) {
                //item_t cur_item;
                if (field.size() == 2)
                    try {
                        // ładujemy obiekt, porównując z haszem md5
                        //cur_item = serialload<item_t > (fs::path(base_dir/field[0]), field[1]);
                        //results.push_back(cur_item);
                        results.push_back(serialload<item_t > (fs::path(base_dir/field[0]), field[1]));
                    } catch (exception::file_not_found e) {
                        sqlite_log << "get(): object file not found: " << field[0] << std::endl;
                    } catch (exception::wrong_md5 e) {
                        sqlite_log << "get(): md5 mismatch: " << field[1] << std::endl;
                    }
            }
            sqlite_log << "get(): restored " << results.size() << " objects\n";
            return results;
        };

    public:

        /**
         * @param dbfile plik bazy danych
         * @param dir bazowy katalog bazy danych
         */
        base(const fs::path & dbfile, const fs::path & dir) :
        base_dir(dir) {
            int result = 0;
	    if(!fs::exists(dbfile.parent_path()))
		fs::create_directories(dbfile.parent_path());
            if(!fs::exists(base_dir))
                fs::create_directories(base_dir);


            result = sqlite3_open(dbfile.string().c_str(), &sqlite_db);
            if (result != SQLITE_OK) {
                sqlite3_close(sqlite_db);
                sqlite_log << sqlite3_errmsg(sqlite_db) << std::endl;
                //std::cout << sqlite_log.str() << std::endl;
                throw exception::sqlite_failed_opening_db();
            }
            try {
                command("CREATE TABLE IF NOT EXISTS objects (filename TEXT, date DATETIME, md5 TEXT PRIMARY KEY)");
            } catch (exception::sqlite_command_error e) {
                //std::cout << sqlite_log.str() << std::endl;
                throw exception::sqlite_failed_creating_db();
            }

        }

        ~base() {
            sqlite3_close(sqlite_db);
        }

        /**
         * @brief podejmowanie obiektów z bazy na podstawie zapytania
         * @param wheres lista warunków typu WHERE wielkość='wartość'
         * @return wektor obiektów spełniających zapytanie
         *
         * Do generowania wheres używamy funkcji where w składni:
         * @code
         * get<item_t>(where(a,b)(c,d)(e,f));
         * @endcode
         * gdy chcemy wyjąć obiekty dla których a=b,c=d i e=f
         */
        template<class item_t>
        std::vector<item_t> get(const std::vector<pair_t> & wheres) {
            std::stringstream s;
            s << "SELECT filename,md5 from objects\nWHERE ";
            foreach(pair_t p, wheres) {
                std::vector<std::string> typeandvalue;
                boost::split(typeandvalue,p.second,boost::is_any_of(":"));

                s << p.first << "=\'" << typeandvalue[1] << "\'";
                if (p != wheres.back())
                    s << " AND ";
            }
            std::vector<item_t> results;
            #pragma omp critical
            results = eval_get<item_t>(s.str());
            return results;
        }

        /**
         * @brief podejmowanie obiektów z bazy na podstawie zapytania
         * @param betweens lista warunków typu WHERE a BETWEEN b AND c
         * @return wektor obiektów spełniających zapytanie
         *
         * Do generowania betweens używamy funkcji between w składni:
         * @code
         * get<item_t>(between(a,b,c)(d,e,f)
         * @endcode
         * gdy chcemy wyjąć obiekty, dla których a należy do [b,c) i
         * d należy do [e,f).
         */
        template<class item_t>
        std::vector<item_t> get(const std::vector<tween_t> & betweens){
            std::stringstream s;
            s << "SELECT filename,md5 from objects\nWHERE ";
            foreach(tween_t p, betweens){
                s << p.field << " BETWEEN \'" << p.begin << "\' AND \'" << p.end << "\'";
                if(p != betweens.back())
                    s << " AND ";
            }
            std::vector<item_t> results;
            #pragma omp critical
            results = eval_get<item_t>(s.str());
            return results;
        }

        /**
         * @brief podejmowanie obiektów z bazy na podstawie zapytania, wersja hybrydowa
         * @param wheres lista warunków typu WHERE wielkość='wartość'
         * @param betweens lista warunków typu WHERE a BETWEEN b AND c
         * @return wektor obiektów spełniających zapytanie
         *
         * Funkcja działa jak logiczna suma funkcji jednoargumentowych get.
         */
        template<class item_t>
        std::vector<item_t> get(const std::vector<pair_t> & wheres, const std::vector<tween_t> & betweens){
            std::stringstream s;
            s << "SELECT filename,md5 from objects\nWHERE ";
            foreach(pair_t p, wheres) {
                std::vector<std::string> typeandvalue;
                boost::split(typeandvalue,p.second,boost::is_any_of(":"));

                s << p.first << "=\'" << typeandvalue[1] << "\'";
                if (p != wheres.back())
                    s << " AND ";
            }
            if(betweens.size()!=0) s << " AND ";
            foreach(tween_t p, betweens){
                s << p.field << " BETWEEN \'" << p.begin << "\' AND \'" << p.end << "\'";
                if(p != betweens.back())
                    s << " AND ";
            }
            std::vector<item_t> results;
            #pragma omp critical
            results = eval_get<item_t>(s.str());
            return results;
        }

        /**
         * \brief zapisuje obiekt w bazie danych SQLite
         *
         * @param item          obiekt do zapisania
         * @param column_values szczegółowe metadane na temat obiektu, wektor par (kolumna,wartość)
         * @param base_dir          ścieżka do katalogu bazowego
         */
        template<class item_t>
        void store(item_t & item, const std::vector<pair_t> & column_values) {
            std::stringstream s1, s2; // dwie części polecenia SQL, w jednej nazwy kolumn, w drugiej wartości
            std::stringstream f; // nazwa pliku
            std::string md5 = md5gen<item_t > (item); // unikalny hash obiektu
            pt::ptime local_time = pt::microsec_clock::local_time(); // aktualny czas z dokładnością do us
            //f << base_dir << "/";

            // tworzymy zapytanie
            s1 << "INSERT INTO objects ( ";
            s2 << " VALUES ( ";

            // sortujemy kolumny alfabetycznie
            std::vector<pair_t> sorted_column_values = column_values;
            std::sort(boost::begin(sorted_column_values), boost::begin(sorted_column_values), pair_t_pred);

            foreach(pair_t & p, sorted_column_values) {
                // uzupełniamy niesitniejące kolumny
                std::vector<std::string>    typeandvalue;
                boost::split(typeandvalue,p.second,boost::is_any_of(":"));

                std::stringstream create;
                create << "ALTER TABLE objects ADD COLUMN " << p.first << " ";
                if(typeandvalue[0]==typeid(std::string).name())
                    create << "TEXT";
                if(typeandvalue[0]==typeid(double).name())
                    create << "REAL";
                if(typeandvalue[0]==typeid(int).name())
                    create << "REAL";
                if(typeandvalue[0]==typeid(float).name())
                    create << "REAL";
                try {
                    command(create.str());
                } catch (exception::sqlite_command_error e) {
                } // SQLite będzie wywalać błąd, dlatego go ignorujemy

                f << p.first << "=" << typeandvalue[1] << ",";
                s1 << p.first << ',';
                s2 << '\'' << typeandvalue[1] << "\' ,";
            }
            s1 << "filename, date, md5 ) ";
            f << "md5=" << md5 << ".z";
            s2 << '\'' << f.str() << "\' , \'" << pt::to_simple_string(local_time) << "\', \'" << md5 << "\' )";

            // wykonanie polecenia SQL
            #pragma omp critical
            command((s1.str() + s2.str()));
            // zapisanie obiektu do pliku
            serialsave<item_t > (item, fs::path(base_dir/f.str()));
        }

        ///dostęp do loga SQLite
        const std::stringstream & log() const {
            return sqlite_log;
        }
    };
};



#endif	/* _BB_SQLITE_H */

