#include "bandit/bandit.h"
#include "Lattice.h"
#include "std.h"

using namespace bandit;

go_bandit([]()
{



    describe("Lattice", []()
    {
        Settings set;
        set.lattice.L = 5;
        set.lattice.W = 5;
        set.lattice.H = 5;
        Lattice l(set);

        it("should copy", [&]()
        {
            Lattice l2 = l;
            //AssertThat(l2, Equals(l));
        });

        it("should parallel copy", [&]()
        {
            std::vector<std::thread> threads;
            for(int i = 0; i < 20; i++)
                threads.push_back(std::thread([&]
                {
                    Lattice l2 = l;
                }));
            for(std::thread & t : threads)
                t.join();
        });

    });

});


int main(int argc, char* argv[])
{
    return bandit::run(argc, argv);
}

