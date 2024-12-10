CC = mpicc
CFLAGS = -std=c99
TARGET = dbscan
SRC = dbscan.c
OUTPUT = results.txt
PROCESSES = 1 2 4 8 16 32 64
RUNS = 3

.PHONY: all clean run average

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) -lm

run: $(TARGET)
	@rm -f $(OUTPUT)
	@echo "Running tests..." >> $(OUTPUT)
	@for np in $(PROCESSES); do \
		echo "Processes: $$np" >> $(OUTPUT); \
		total_time=0; \
		for i in `seq 1 $(RUNS)`; do \
			echo "Run $$i with $$np processes..."; \
			time_output=$$(mpirun -np $$np ./$(TARGET) | grep "Time for clustering" | awk '{print $$4}'); \
			echo "Run $$i: $$time_output seconds" >> $(OUTPUT); \
			total_time=$$(echo $$total_time + $$time_output | bc); \
		done; \
		average_time=$$(echo "scale=6; $$total_time / $(RUNS)" | bc); \
		echo "Average time for $$np processes: $$average_time seconds" >> $(OUTPUT); \
		echo "==============================" >> $(OUTPUT); \
	done

average: run
	@echo "Results saved in $(OUTPUT)"

clean:
	rm -f $(TARGET) $(OUTPUT)
