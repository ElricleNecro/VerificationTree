########################################################
#	Variables de compilations :
########################################################

#!|-----------------------------------------------------|
#!|			Variables			|
#!|-----------------------------------------------------|

#!CC :
#!	Compilateur à utiliser.
#!
CC=gcc

CFLAG=-std=c99 -O3 -W -Wall -Wshadow -Wcast-qual \
      -Wcast-align -Wsign-compare -Wstrict-prototypes \
      -Wredundant-decls \
      -Wnested-externs \
      -ffloat-store -Wunreachable-code -Wwrite-strings
#-ggdb -Wmissing-declarations
EXTRA=-pg
#-pg

#!IOPERSO :
#!	IOPERSO = 1 :Utilisation des fonctions personnelles de lecture des fichiers Gadget, plutôt que d'utiliser celles de Springel améliorer.
#!	Défaut : IOPERSO = ""
#!

#!DEBUG :
#!	Flag de Debug à passer au compilateur.
#!

#!TIMER :
#!	TIMER = -DUSE_TIMER : Affiche le temps d'éxecution des fonctions CalcPotentiel et DensityCenter.
#!	TIMER = ""          : N'affiche aucun timer.
#!	Valeur par défaut   : -DUSE_TIMER.
#!
TIMER=-DUSE_TIMER

DEBUG+=-D__DEBUG_VOIS_LOG
DEBUG+=-DUSE_STRUCT_PART
DEBUG+=-DTREE_CM_BUILD
DEBUG+=-DUSE_STRUCT_PART

ifeq ($(IOPERSO),1)
DEBUG+=-DIO_PERSO
endif

#-DP_DBG_TREECODE_P_CALC2 -DTREE_CM_BUILD_DEBUG_ -DTREE_CALCPOT_DEBUG_ -DPREV_NB_WEAK_P -DPREV_NB_NEG_P -DPREV_LISSAGE_P -DUSE_MCHECK
# -DP_DBG_TREECODE_P_CALC -DP_DBG_TREECODE_P
#-DP_DBG_TREECODE_P_CALC
#-DP_DBG_TREECODE_P -DP_DBG_TREECODE_P_CALC

INC=-I $$HOME/.C_C++/include -I $$HOME/.local/include
LFLAG=-L $$HOME/.C_C++/lib -I $$HOME/.local/lib
LINK=-lm

#!PREFIX :
#!	Dossier dans lequel sera installé le programme.
#!
PREFIX=$$HOME/.local/

########################################################
#	Executable divers :
########################################################
RM=/bin/rm
VAL=/usr/bin/valgrind
GDB=/usr/bin/gdb

########################################################
#	Variables d'archivage :
########################################################
TAR=tar -cvzf
CTAR=VerifTreeCode
EXTT=tar.gz

########################################################
#	Sources :
########################################################
MAIN=Verif.c
MAIN2=BruteForce.c

ifeq ($(IOPERSO),1)
SRC=Verif_tools.c tree.c iogadget2.c utils.c types.c
else
SRC=Verif_tools.c tree.c utils.c types.c io_snapshot.c
endif

SRC_TREE=#tree_create.c tree_voisin.c
#commun.c Voisin1.c

OBJ=$(SRC:.c=.o) $(SRC_TREE:.c=.o)
HEADERS=$(SRC:.c=.h)

########################################################
#	Executable :
########################################################
#!EXEC :
#!	Nom de l'éxecutable principale (défaut : Verif).
#!
EXEC=$(MAIN:.c=)

#!EXEC2 :
#!	Nom de l'éxecutable brute-force.
#!
EXEC2=$(MAIN2:.c=)

#!|-----------------------------------------------------|
#!|			Cibles				|
#!|-----------------------------------------------------|

.PHONY : clean clean-all doc install tar gdb valgrind

#!all-single :
#!	Crée l'éxecutable du programme.
#!
all-single:$(EXEC)

#!all :
#!	Crée à la fois l'éxecutable du projet et la version brute-force du code.
#!
all:$(EXEC) $(EXEC2)

$(MAIN:.c=.o):$(MAIN)
	$(CC) $(CFLAG) $(EXTRA) $(DEBUG) $(DBGFLAG) $(TIMER) $(INC) -c $<

$(MAIN2:.c=.o):$(MAIN2)
	$(CC) $(CFLAG) $(EXTRA) $(DEBUG) $(DBGFLAG) $(TIMER) $(INC) -c $<

tree_create.o:tree_create.c tree.h
	$(CC) $(CFLAG) $(EXTRA) $(DEBUG) $(DBGFLAG) $(TIMER) $(INC) -c $<

tree_voisin.o:tree_voisin.c tree.h
	$(CC) $(CFLAG) $(EXTRA) $(DEBUG) $(DBGFLAG) $(TIMER) $(INC) -c $<

%.o:%.c %.h
	$(CC) $(CFLAG) $(EXTRA) $(DEBUG) $(DBGFLAG) $(TIMER) $(INC) -c $<

$(EXEC): $(MAIN:.c=.o) $(OBJ) $(HEADERS)
	$(CC) $(CFLAG) $(EXTRA) $(DEBUG) $(DBGFLAG) $(INC) $(LFLAG) $(OBJ) $(MAIN:.c=.o) -o $@ $(LINK)

$(EXEC2): $(MAIN2:.c=.o) $(OBJ) $(HEADERS)
	$(CC) $(CFLAG) $(EXTRA) $(DEBUG) $(DBGFLAG) $(INC) $(LFLAG) $(OBJ) $(MAIN2:.c=.o) -o $@ $(LINK)

test_temp:
	$(CC)  $(CFLAG) $(EXTRA) $(DEBUG) $(DBGFLAG) $(TIMER) ../Tree_Code/rand.c utils.c tree.c Verif_tools.c types.c test_Temp.c -o $(LFLAG) test_temp  $(LINK)

#!help :
#!	Affiche cette aide.
#!
help:
	@grep "^#!" Makefile | sed -e 's:^#!::'

########################################################
#	Debug :
########################################################
#!valgrind :
#!	Lance l'éxecutable principal dans valgrind pour vérifier la mémoire, en utilisant les options données dans la variables ARGS.
#!
valgrind: $(EXEC)
	$(VAL) --tool=memcheck --num-callers=40 --leak-check=full --track-origins=yes --show-reachable=yes ./$(EXEC) $(ARGS)

#!gdb :
#!	Lance l'éxecutable dans gdb pour déboggage.
#!
gdb: $(EXEC)
	$(GDB) $(EXEC)

########################################################
#	Installation :
########################################################
#!install :
#!	Installe le programme dans $PREFIX/bin.
#!
install:$(EXEC)
	@mkdir -p $(PREFIX)/bin && cp $(EXEC) $(PREFIX)/bin/.

########################################################
#	Archivage :
########################################################
#!tar :
#!	Compresse les sources dans une archive tar.gz.
#!
tar:
	@$(TAR) $(CTAR).$(EXTT) $(SRC) $(MAIN) $(HEADERS) Makefile cte_phys.h

########################################################
#	Documentation :
########################################################
#!doc :
#!	Génére la doc des fonctions.
doc:
	doxygen Doxyfile

########################################################
#	Nettoyage :
########################################################
#!clean :
#!	Nettoie les fichiers objets et back-up.
#!
clean:
	$(RM) -rf *.o *.c~

#!clean-all :
#!	Nettoie tout : fichiers objets, back-up, documentation, éxecutable, ...
#!
clean-all:clean
	$(RM) -rf $(EXEC) Doc/

