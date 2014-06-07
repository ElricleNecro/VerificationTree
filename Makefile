########################################################
#	Variables de compilations :
########################################################

#!|-----------------------------------------------------|
#!|		Variables de compilations		|
#!|-----------------------------------------------------|

#!CC :
#!	Compilateur à utiliser.
#!
CC=gcc

#!|-----------------------------------------------------|
#!|		   Variables de debug			|
#!|-----------------------------------------------------|

#!TIMER :
#!	TIMER = -DUSE_TIMER : Affiche le temps d'éxecution des fonctions CalcPotentiel et DensityCenter.
#!	TIMER = ""          : N'affiche aucun timer.
#!	Valeur par défaut   : -DUSE_TIMER.
#!
TIMER=-DUSE_TIMER

#!SQLITE3 :
#!	Utilisation de base de donnée plutôt que des fichiers.
#!	Valeur par défaut : -DUSE_SQLITE3
#SQLITE3=-DUSE_SQLITE3 $(shell pkg-config --cflags sqlite3)

#!HDF5 :
#!	Utilisation du système de fichier HDF5 pour l'enregistrement
#!	Valeur par défaut : -DUSE_HDF5 -lhdf5
HDF5=-DUSE_HDF5

#!DEBUG :
#!	Flag de Debug à passer au compilateur.
#!
#!		(x) __DEBUG_VOIS_LOG                  : Pour le calcul du centre de densité : affiche la densité locale autour de chaque particule.
#DEBUG+=-D__DEBUG_VOIS_LOG
#!		(x) USE_STRUCT_PART                   : Active l'utilisation d'une structure à la Gadget plutôt que d'un tableau 2D pour contenir les infos sur les particules. DÉPRETIÉ.
DEBUG+=-DUSE_STRUCT_PART
#!		(X) TREE_CM_BUILD                     : Active le calcul du centre de gravité en même temps que l'arbre se construit.
#DEBUG+=-DTREE_CM_BUILD
#!		(x) TEST_VOISIN_LOG_SHUFFLE           : Calcul une première fois le centre, re-mélange les particules et recalcul à nouveau. Si les résultats sont les mêmes, tout est bon, sinon l'arbre ou
#!					      		l'algorithme de recherche des voisins ne fonctionne pas correctement.
#DEBUG+=-DTEST_VOISIN_LOG_SHUFFLE
#!		(x) TEST_INFLUENCE_MODIF_MARCHE_ARBRE : Fait deux fois de suite le calcul du centre afin de voir l'influence des itérations du calcul. Les différents résultats devraient à peine bouger.
#DEBUG+=-DTEST_INFLUENCE_MODIF_MARCHE_ARBRE
#!		(x) PERIODIC                          : active la prise en compte des conditions périodiques lors du centrage de l'objet.
#DEBUG+=-DPERIODIC
#!		(x) DOUBLE_BOUCLE                     : active le calcul des voisins en utilisant une boucle remplissant un premier tableau classé par distance à la particule test puis une boucle insérant
#!				    			les particules de ce tableau dans le tableau des voisins.
#DEBUG+=-DDOUBLE_BOUCLE
#!		(x) USE_VOIS_QSORT                    : Utilise un qsort pour garder les tableaux trié lors de la recherche des voisins. Ne signifie quelque chose que si DOUBLE_BOUCLE est activé.
#DEBUG+=-DUSE_VOIS_QSORT
#-DP_DBG_TREECODE_P_CALC
#!		(x) TREE_CALCPOT_DEBUG_               : Affiche des infos de déboguage du potentiel.
#DEBUG+=-DTREE_CALCPOT_DEBUG_
#!		(x) USE_NEWDISTCALC		      : Utilise : fmax( 0., fabs( root->x - part->x ) - root->cote/2.0) pour le calcule de la distance particule--cube.
#DEBUG+=-DUSE_NEWDISTCALC
#DEBUG+=-DUSE_FILE
#!		(x) USE OLDWAY			      : permet de donner directement le type de particule à charger plutôt qu'une puissance de 2 du type (2^type).
#DEBUG+=-DOLDWAY
#!		(x) DBG_NEWWAY			      : quelques affichages pour déboguer le chargement des particules utilisant les bits.
#DEBUG+=-DDBG_NEWWAY

#DEBUG+=-DACTIVATE_FoF
# DEBUG+=-DACTIVATE_SPHERICAL_SELECTION
# DEBUG+=-DACTIVATE_SPHERICAL_SELECTION_AFTER
#DEBUG+=-DDEBUG_FOF
DEBUG+=-D__DENSITYCENTER_NOPROGRESS_P
DEBUG+=-D__POTENTIALCENTER_NOPROGRESS_P
DEBUG+=-D__FoF_NOPROGRESSBAR_P

#DEBUG+=-DUSE_MCHECK

#!IOPERSO :
#!	IOPERSO = 1 :Utilisation des fonctions personnelles de lecture des fichiers Gadget, plutôt que d'utiliser celles de Springel améliorer.
#!	Défaut : IOPERSO = ""
#!
ifeq ($(IOPERSO),1)
DEBUG+=-DIO_PERSO
endif

#-DP_DBG_TREECODE_P_CALC2 -DTREE_CM_BUILD_DEBUG_ -DTREE_CALCPOT_DEBUG_ -DPREV_NB_WEAK_P -DPREV_NB_NEG_P -DPREV_LISSAGE_P -DUSE_MCHECK
# -DP_DBG_TREECODE_P_CALC -DP_DBG_TREECODE_P
#-DP_DBG_TREECODE_P_CALC
#-DP_DBG_TREECODE_P -DP_DBG_TREECODE_P_CALC

EXTRA=
#-pg
INC=-I $$HOME/.local/include -I include/ -I octree/ -I /softs/hdf5/1.8.11/include

CFLAG+=-std=c99 -O3 -W -Wall -Wshadow -Wcast-qual \
       -Wcast-align -Wsign-compare -Wstrict-prototypes \
       -Wredundant-decls \
       -Wnested-externs \
       -ffloat-store -Wunreachable-code -Wwrite-strings \
      $(DEBUG) $(DBGFLAG) $(TIMER) $(SQLITE3) $(HDF5) $(EXTRA) $(INC) -g3
#-ggdb -Wmissing-declarations
#-floop-nest-optimize

ifeq ($(shell hostname),"Archlinux-Dell")
	CFLAG+=-fsanitize=address
endif

LINK=-lm $(shell pkg-config --libs sqlite3) -lhdf5 -loctree
LFLAG=-L $$HOME/.local/lib -L lib/ -L /softs/hdf5/1.8.11/lib

#!|-----------------------------------------------------|
#!|		  Variables de dossiers			|
#!|-----------------------------------------------------|

#!PREFIX :
#!	Dossier dans lequel sera installé le programme.
#!
PREFIX=$$HOME/.local/

SRCDIR=src
INCDIR=include
OBJDIR=build
EXECDIR=build

#!|-----------------------------------------------------|
#!|		Variables d'Éxecutable			|
#!|-----------------------------------------------------|

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
SRC=Verif_tools.c tree.c utils.c types.c io_snapshot.c HDF5.c
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
EXEC=$(EXECDIR)/$(MAIN:.c=)

#!EXEC2 :
#!	Nom de l'éxecutable brute-force.
#!
EXEC2=$(EXECDIR)/$(MAIN2:.c=)

LIB_OCT=liboctree.so

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
all: $(EXEC) $(EXEC2)

$(OBJDIR)/$(MAIN:.c=.o):$(SRCDIR)/$(MAIN) Makefile lib/$(LIB_OCT)
	@echo -e "\033[32m----------------------------------------------------------\033[00m"
	@echo -e "\033[31mCompiling " $< "\033[00m"
	$(CC) $(CFLAG) -c $< -o $@

$(OBJDIR)/$(MAIN2:.c=.o):$(SRCDIR)/$(MAIN2) Makefile lib/$(LIB_OCT)
	@echo -e "\033[32m----------------------------------------------------------\033[00m"
	@echo -e "\033[31mCompiling " $< "\033[00m"
	$(CC) $(CFLAG) -c $< -o $@

tree_create.o:tree_create.c tree.h
	@echo -e "\033[32m----------------------------------------------------------\033[00m"
	@echo -e "\033[31mCompiling " $< "\033[00m"
	$(CC) $(CFLAG) -c $<

tree_voisin.o:tree_voisin.c tree.h
	@echo -e "\033[32m----------------------------------------------------------\033[00m"
	@echo -e "\033[31mCompiling " $< "\033[00m"
	$(CC) $(CFLAG) -c $<

$(OBJDIR)/%.o:$(SRCDIR)/%.c $(INCDIR)/$(HEA) Makefile
	@echo -e "\033[32m----------------------------------------------------------\033[00m"
	@echo -e "\033[31mCompiling " $< "\033[00m"
	$(CC) $(CFLAG) -c $< -o $@

$(EXEC):$(OBJDIR)/$(MAIN:.c=.o) $(foreach x, $(OBJ), $(OBJDIR)/$(x))
	@echo -e "\033[32m----------------------------------------------------------\033[00m"
	@echo -e "\033[31mBuilding " $@ "\033[00m"
	$(CC) $(CFLAG) $(LFLAG) $^ -o $@ $(LINK)

$(EXEC2):$(OBJDIR)/$(MAIN2:.c=.o) $(foreach x, $(OBJ), $(OBJDIR)/$(x))
	@echo -e "\033[32m----------------------------------------------------------\033[00m"
	@echo -e "\033[31mBuilding " $@ "\033[00m"
	$(CC) $(CFLAG) $(LFLAG) $^ -o $@ $(LINK)

lib/$(LIB_OCT):$(OBJDIR)/octree.o
	@mkdir -p lib
	${CC} -shared $(LFLAG) $< -o $@

$(OBJDIR)/octree.o:octree/octree.c
	$(CC) $(CFLAG) -fPIC -c $< -o $@

test_temp:
	$(CC)  $(CFLAG) ../Tree_Code/rand.c utils.c tree.c Verif_tools.c types.c test_Temp.c -o test_temp $(LFLAG)

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
	@mkdir -p $(PREFIX)/bin
	@install -m 755 $(EXEC) $(PREFIX)/bin
	@mkdir -p $(PREFIX)/lib
	@install -m 755 lib/$(LIB_OCT) $(PREFIX)/lib
#@mkdir -p $(PREFIX)/bin && cp $(EXEC) $(PREFIX)/bin/.

########################################################
#	Archivage :
########################################################
#!tar :
#!	Compresse les sources dans une archive tar.gz.
#!
tar:
	@$(TAR) $(CTAR).$(EXTT) $(foreach x, $(SRC), $(SRCDIR)/$(x)) $(SRCDIR)/$(MAIN) $(foreach x, $(HEADERS), $(INCDIR)/$(x)) Makefile $(INCDIR)/cte_phys.h octree

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
	$(RM) -rf $(OBJDIR)/*.o $(SRCDIR)/*.c~

#!clean-all :
#!	Nettoie tout : fichiers objets, back-up, documentation, éxecutable, ...
#!
clean-all:clean
	$(RM) -rf $(EXEC) Doc/

