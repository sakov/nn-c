all:
	./test.sh
mpi:
	./mpitest.sh
cmp:
	@echo -n "  running mpitest.sh..."
	@time -f " %es" sh -c '"$$0" >/dev/null 2>&1' ./mpitest.sh
	@for file in lin.txt nn-inf.txt nn-0.txt nn-ns.txt; do mv $$file $${file}.mpi; done
	@echo -n "  running test.sh..."
	@time -f " %es" sh -c '"$$0" >/dev/null 2>&1' ./test.sh
	@echo "  comparing results:"
	@for file in lin.txt nn-inf.txt nn-0.txt nn-ns.txt; do echo -n "    "; diff -qs $$file $${file}.mpi; done
clean:
	rm -f lin*.txt* nn-*.txt* *~ core
