CC = gcc
CFLAGS = -O3
COMMON_HDR = defines.h
COMMON = compute.o read.o scatdb.o print.o /usr/lib64/libm.so

all : xrd_calc xrd_calc_alt test_main find_osf search_hbs_f search_f canvas_hbs_f merge_same_refs

canvas_hbs_f : canvas_hbs_f.o $(COMMON)

search_f : search_f.o $(COMMON)

search_hbs_f : search_hbs_f.o $(COMMON)

xrd_calc : xrd_calc.o $(COMMON)

xrd_calc_alt : xrd_calc_alt.o $(COMMON)

test_main : test_main.o $(COMMON)

find_osf : find_osf.o $(COMMON)

merge_same_refs : merge_same_refs.o $(COMMON)

xrd_calc.o : $(COMMON_HDR)

xrd_calc_alt.o : $(COMMON_HDR)

test_main.o : $(COMMON_HDR)

find_osf.o : $(COMMON_HDR)

search_f.o : $(COMMON_HDR)

search_hbs_f.o : $(COMMON_HDR)

canvas_hbs_f.o : $(COMMON_HDR)

merge_same_refs.o : $(COMMON_HDR)

$(COMMON) : $(COMMON_HDR)

clean :
	rm merge_same_refs xrd_calc xrd_calc_alt test_main find_osf search_f search_hbs_f canvas_hbs_f merge_same_refs.o xrd_calc.o xrd_calc_alt.o test_main.o find_osf.o search_f.o search_hbs_f.o canvas_hbs_f.o compute.o read.o scatdb.o print.o

