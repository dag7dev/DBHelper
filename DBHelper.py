# dag7 - 2021 - GPLv3

import json
import math
import os
import re
import sys
import copy
import getopt

verbose = False      # stampa qualcosa in più per capire il procedimento e la logica
very_verbose = False # stampa davvero tutto

# flag
isam_opt = False
btree_opt = False
hash_opt = False

# cartelle
jsondir = "json"
phyjsondir = "phyjson"

########
# UTILS
########
# crea una stringa per stampare in un bel modo un dizionario
def dict_to_str(dictionary):
    s = str()

    for key in dictionary:
        for value in dictionary[key]:
            s += key + " -> " + value + ", "

    # rimuoviamo la virgola alla fine
    return s[:-2]

# stampa l'header con tutti i dati
def print_header(R, F, keys, ro, X):
    s = "F:\t\t"

    print("Launched with parameters:")
    print("R:\t\t" + str("".join(sorted(list(R)))))

    print("F:\t\t" + dict_to_str(F))
    print("Keys (if any):\t" + str(keys) if len(keys) != 0 else "Keys (if any):\t--no key given--")
    print("Ro (if any):\t" + str(ro) if len(ro) != 0 else "Ro (if any):\t--no ro given--")
    print("X (if any):\t" + str(X) if len(X) != 0 else "X (if any):\t--no X given--")
    print("\n")

# stampa un errore
def print_err(arr, var_name):
    if not arr:
        input("! WARNING: you are testing an empty " + var_name + ". Set up " + var_name + " in your json for best results.\nPress ENTER to continue...")

# stampa un messaggio di help
def usage(extra=False):
    print("Syntax: DBHelper.py <option> <parameter:json_file_containing_a_scheme | all>")
    print("<option>:")
    print("\t-h --help: \n\t\tdisplays this message.")
    print("\t-v --verbose: \n\t\tshows `what's going on` messages")
    print("\t-x --extra-verbose: \n\t\tshows anything that this program could print")
    print("\t-n --third-nf: \n\t\ttells if a scheme is in 3NF or not")
    print("\t-k --find-keys: \n\t\tfind and shows keys of a scheme")
    print("\t-t --test-keys: \n\t\ttest keys defined in the scheme, if they are keys or not")
    print("\t-c --closure: \n\t\tclosure of a given X attributes")
    print("\t-m --minimal: \n\t\tfind a minimal coverage of a scheme")
    print("\t[WIP] -p --preserve: \n\t\tcheck if decompositions in ro preserve FDs --> currently it states the FDs you manually need to check")
    print("\t[WIP] -l --lossless: \n\t\tcheck if R, given decompositions in ro, has a lossless join. --> currently it prints the big table at its initial state.")
    print("\t--json: \n\t\tcustom json folder (default: .)")
    
    if extra:
        print("\nPHYSICAL DATA ORGANIZATION")
        print("\t--physical-isam: \n\t\tsolve a json containing isam exercise")
        print("\t--physical-hash: \n\t\tsolve a json containing hash exercise")
        print("\t[WIP]--physical-btree: \n\t\tsolve a json containing btree exercise")
        print("\t--phyjson: \n\t\tcustom json folder (default: phyjson)")
    else:
        print("More options are available, type `python DBHelper.py --extended-help` for more")
    print("NOTE: WIP functions are WIP for a reason: you're welcome to give a hand where needed!")
    
    exit(1)

########
# JSON
########
# controlla se è un json valido
def validate_json(json_data):
    try:
        json.loads(json_data)
    except ValueError:
        return False
    return True

# array may contains a string or an array, depending on user
def load_arr_or_str(key, arr):
    final = set()

    # check if exists and is not null
    if key in arr and key != "":
        if type(arr[key]) == str:   # se unica stringa
            final = set(re.split(', ', arr[key]))
        else:
            if len(arr[key]) == 1:
                final = set(re.split(', ', arr[key][0]))
            else:
                for r in arr[key]:    # altrimenti l'utente si è sbagliato ed è una lista
                    final.add(r)
    else:
        final = {}
    return final

# parsa i dati da un json di un esercizio
def parse_json_data(json_data):
    '''
        How JSON is made:
        R: string -> schema, it contains a group of attributes
        F: list containing dicts -> each key: value contains a list of associated dependencies
        keys: optional -> contains a list of keys (made up by strings)
        ro: optional -> contains possible decompositions 
        X: optional -> contains subsets of R to test
    '''
    R = set()
    F = dict()
    X = list()
    keys = list()
    ro = list()
    
    # check if F exists
    if not json_data["F"]:
        print("Couldn't retrieve F from json. Verify that F is present and try again")
        exit(1)

    # get F+: the functional dependencies are processed AND decomposed "BY DEFAULT"
    for df in json_data["F"]:
        splitted_df = re.split(' -> ', df)
        
        left_side = splitted_df[0]
        right_side = splitted_df[1]

        if len(left_side) != 1:
            if left_side not in F:  # casi come A -> BC
                F[left_side] = list(right_side)
            else:
                for r in right_side:
                    if r not in F[left_side]:
                        F[left_side].append(r)
        else:   # casi come A -> BC
            if left_side not in F:
                F[left_side] = list(right_side[0])
                sliced = right_side[1:]
                # casi come AB -> CD
                for r in sliced:
                    F[left_side].append(r)
            else:
                for r in right_side:
                    F[left_side].append(r)

     # check if R exists
    if not json_data["R"]:
        print("R has not been provided. It will be calculated!")
        left = set()
        right = set()
        
        for f in F.keys():
            left = left.union(set(f))
        
        for r in F.values():
            right = right.union(set(r))
        
        R = left.union(right)
    else:
        # get the schema
        R = set(json_data["R"])

    # retrieval X, ro, keys
    X = load_arr_or_str('X', json_data)
    ro = load_arr_or_str('ro', json_data)
    keys = load_arr_or_str('keys', json_data)
    
    return R, F, sorted(list(set(keys))), set(sorted(list(ro))), set(sorted(list(X))) # converting to set and then to list, to avoid duplicates

# parsa i dati da un json
def parse_physical_json(json_data):
    # isam
    if json_data['type'].upper() == "ISAM" and isam_opt:
        if "NR" not in json_data or "R" not in json_data or "K" not in json_data or "CB" not in json_data or "P" not in json_data or not json_data["NR"] or not json_data["R"] or not json_data["K"] or not json_data["CB"] or not json_data["P"]:
            print("Missing data from json. Verify that data is present and it is correct and try again")
            exit(1)

        return json_data['type'], json_data["NR"], json_data["R"], json_data["K"], json_data["CB"], json_data["P"], json_data["filled_space"] if json_data["filled_space"] else 0, False if not json_data["floor_min"] else True
    
    # hash
    elif json_data['type'].upper() == "HASH" and hash_opt:
        if "NR" not in json_data or "R" not in json_data or "CB" not in json_data or "P" not in json_data or "B" not in json_data or not json_data["NR"] or not json_data["R"] or not json_data["CB"] or not json_data["P"] or not json_data["B"]:
            print("Missing data from json. Verify that data is present and try again")
            exit(1)
        return json_data['type'], json_data["NR"], json_data["R"], 0 if "K" not in json_data or not json_data["K"] else json_data["K"] , json_data["CB"], json_data["P"], 0 if "filled_space" not in json_data or not json_data["filled_space"] else json_data["filled_space"], False if "floor_min" not in json_data or not json_data["floor_min"] else True, json_data["B"], 0 if "no_access" not in json_data or not json_data["no_access"] else json_data["no_access"]
    
    elif json_data['type'].upper() == "BTREE" and btree_opt:
        print("COMING SOON!")
        exit(1)

    else:
        print("Filetype seems to be: " + json_data['type'])
        print("Type hasn't been specified in the json or it is wrong <ISAM, BUCKET, HASH, BTREE>")
        print("Make sure to correctly specify the type of test and try again.")
        exit(1)
        
# carica un json
def json_load(jsondir, json_path, physical = False):
    global verbose
    global very_verbose

    if not jsondir.endswith("/"):
        jsondir = jsondir + '/'

    json_path = jsondir + json_path

    # safe parse filename: append .json if missing (a user could miss json extension)
    if(not json_path.endswith(".json")):
        json_path = json_path + ".json"

    print("File:\t\t" + str(json_path))

    # natural domain
    if(not os.path.exists(json_path) or     # check if file exists 
        not os.path.isfile(json_path)):       # check if it's a file
        print(json_path + " file doesn't exist or it is not a file. Make sure it is a regular json and try again.")
        exit(1)

    # load file
    try:
        with open(json_path) as f:
            data = json.load(f)
            if not physical:
                return parse_json_data(data)
            else:
                return parse_physical_json(data) 
    except ValueError as ve:
        print(ve)
        print(json_path + " it is not a regular json! Build it up according to the specs and try again")

###########
# MODALITA
###########
# utility: calcola chiusura di X+F 
def closure(R, F, X):
    '''
        Z = X
        S = { X -> Y app F ^ Y è a destra ^ Y è contenuto in Z }
        while non ho aggiunto qualcosa di nuovo:
            Z = Z U S
            S = { X -> Y app F ^ Y è a destra ^ Y è contenuto in Z }

        L'algoritmo termina quando non aggiunge qualcosa di nuovo
    '''
    F_copy = copy.deepcopy(F)
    X = set(X)
    
    Z = X
    S = set()
    added_something = False

    # caricamento di S all'inizio
    # scorro il dizionario contenente le df
    for left_side in list(F_copy):
        # se la parte sinistra è un sottoinsieme di X
        # allora vado a vedere le parti destre e me le salvo in S 
        if set(left_side).issubset(set(X)):
            for el in F_copy[left_side]:    # per ogni elemento a destra
                S.add(el)
                added_something = True
            F_copy.pop(left_side)   # cancello da F_Copy, cosi non lo devo riesaminare la prossima volta
    
    # se ho aggiunto qualcosa di nuovo continuo
    while(added_something):
        added_something = False

        Z = Z.union(S)

        # scorro il dizionario contenente le df
        for left_side in list(F_copy):
            # se la parte sinistra è un sottoinsieme di S
            # allora vado a vedere le parti destre e me le salvo in S 
            
            # nel caso di dip composte
            arr_left = []
            for v in left_side:
                arr_left.append(v)

            if set(arr_left).issubset(set(Z)) and len(set(arr_left)) != 0:
                for el in F_copy[left_side]:    # ogni elemento a destra lo aggiungo in F
                    S.add(el)
                    added_something = True
                F_copy.pop(left_side)   # cancello da F_Copy, cosi non lo devo riesaminare la prossima volta
    
    return sorted(list(Z))

# dato un insieme di attributi su F calcola la loro chiusura 
def closures_calc(R, F, X):
    group_of_closures = []
    
    for X in X:
        if very_verbose : print("Detecting closure " + X + "...")
        result = closure(R, F, X)
        group_of_closures.append(result)
        
        if result != []:
            if very_verbose : print("Closure seems to be: " + str(result) + "\n")
        else:
            if very_verbose : print("No closure detected" + "\n")

    return group_of_closures

# calcola la copertura minima
def minimal_coverage(R, F):
    '''
        1. Decomponi
        2. Riduci a sinistra
        3. Riduci a destra

        A sinistra posso averci millecinquecento attributi
        Ogni volta che mi ritrovo roba tipo: AC -> BD

        1. AC -> B, AC -> D
        2. calcolo la chiusura delle robe a sinistra (A e C)
            Verifico se la roba a destra è contenuta in una delle due chiusure della roba a sinistra

            - verifico se B è contenuto in A+F e
            - verifico se B è contenuto in C+F
            A+ = { A, E, D }
            C+ = { C }
            Nessuna delle due chiusure contiene B (roba a destra) e quindi la df va per forza preservata

            - verifico se D è contenuto in A+F
            - verifico se D è contenuto in C+F
            La chiusura di A+ contiene D, e quindi la df diventa A -> D (prima era AC -> D)

        3. a turno elimino...

    '''
    
    # mi ricopio F
    F_simple = copy.deepcopy(F)

    #################################
    # tolgo le ridondanze a sinistra
    #################################
    if verbose : print("Trying to simplify left sides...\n")

    for left_side in copy.deepcopy(F_simple).keys():
        # se a sinistra non ho singleton
        if len(left_side) != 1:
            if verbose : print(left_side + " has been found with multiple attributes in the left side.")
            if very_verbose : print("Calculating closures:")
            closures = closures_calc(R, F, left_side)

            # per ogni df
            right_sides = F_simple[left_side]
            for right_side in right_sides:
                left_side_orig = left_side  # mi salvo la dipendenza originale
                
                i = 0   # contatore delle chiusure
                
                # vai a controllare nella chiusura se la parte sinistra e' nella chiusura
                for c in closures:
                    # se a destra singleton
                    if right_side in c and left_side in F_simple:

                        # rimuovi la dipendenza che sto analizzando
                        if len(right_sides) == 1:
                            # rimuovi la dipendenza che sto analizzando
                            F_simple.pop(left_side)                            
                        else:
                            # se ho un array, rimuovi l'entry se esiste ovviamente
                            if right_side in right_sides:
                                right_sides.pop(right_sides.index(right_side))

                        # e inseriscila semplificata
                        F_simple[left_side_orig[i]].append(right_side)

                    i += 1
    
    # gestisco casi nel quale chiavi contengono altre chiavi con la stessa parte a destra
    # esempio: AB -> D e ABC -> D (C gia' determinato) 
    for key in copy.deepcopy(F_simple).keys():
        for key2 in list((F_simple).keys()): 
            if set(key2).issubset(set(key)) and key != key2:
                # devo cancellarlo se determinano la stessa roba
                i = 0
                for value in F_simple[key]:
                    for value2 in F_simple[key2]:
                        if value == value2:
                            if len(F_simple[key]) == 1:
                                F_simple.pop(key)
                            else:
                                F_simple[key].pop(F_simple[key].index(i))

                    i+=1
                break
        
    if verbose : print("\nSimplified left sides of F: " + dict_to_str(F_simple))

    ##################################
    # tolgo le ridondanze a destra
    ##################################
    # AC -> B proviamo F’= {A -> D, A -> E, B -> E, E -> D }  otteniamo   (AC)+F’={A, C, D,E }
    F_temp = copy.deepcopy(F_simple)

    if verbose : print("\nReducing right sides:")

    for key in F_simple.copy():
        values = F_simple[key]
        for value in values:
            # prima della modifica me lo salvo
            rollback = copy.deepcopy(F_temp)
            if len(F_temp[key]) == 1:
                F_temp.pop(key)
            else:
                F_temp[key].pop(F_temp[key].index(value))
            
            # calcola le chiusure
            cl = closure(R, F_temp, key)
          
            if value in set(cl):
                if verbose : print("Trying to delete " + key + " -> " + value + " OK" + " -> " + value + " is in " + key + "+ = " + str(cl))
                F_simple = copy.deepcopy(F_temp)
            else:
                if verbose : print("Trying to delete " + key + " -> " + value + " NO")
                
                # ripristina la copia prima della modifica
                F_simple = copy.deepcopy(rollback)
                F_temp = copy.deepcopy(rollback)
        
    return F_simple

# funzione: restituisce le chiavi dato R schema e F insieme di dipendenze funzionali
def find_keys(R, F):
    
    min_coverage = minimal_coverage(R, F)
    print("\nMinimal coverage: " + dict_to_str(min_coverage))
    
    # dichiarazione
    keys = []   # array chiavi finali
    keys_f = set()
    values_f = set()

    for key in list(min_coverage.keys()):
        decomposed_set = set(key)
        keys_f = keys_f.union(decomposed_set)
    
    # devo convertire la roba nel set in elementi semplici
    for value in list(min_coverage.values()):
        decomposed_set = set(value)
        values_f = values_f.union(decomposed_set)
    
    
    # essentials, tutti gli attributi che non possono essere derivati da nessuna df
    # common, tutti gli attributi che sono presenti in entrambi i lati di una df
    # useless, solo quelli a destra
    essentials = keys_f.difference(values_f)
    common = keys_f.intersection(values_f)
    useless = values_f.difference(keys_f)

    # in caso non siano contenuti attributi "non determinati", non in F li aggiungo 
    missing_attributes = R.difference(common.union(essentials.union(useless)))
    if missing_attributes != {}:
        for A in missing_attributes:
            essentials.add(A)

    # calcolo la chiusura delle parti
    essentials_cl = set(closure(R, F, list(essentials)))

    # questa è la chiave: determina tutto lo schema
    if essentials_cl == R:
        keys.append("".join(sorted(essentials)))
    else:
        merged_set = set()
        common_sorted = sorted(list(common))
        left_sides = copy.deepcopy(list(min_coverage.keys()))
        for left_side in left_sides:
            merged_set = essentials.union(set(left_side))
            merged_cl = set(closure(R, F, list(merged_set)))
            
            if merged_cl == R:
                keys.append("".join(sorted(merged_set)))
        
    return sorted(list(set(keys)))

# dato K insieme di chiavi, controlla se qualcosa in K è una chiave
def are_keys(R, F, cand_keys):
    str_output = ""

    keys = find_keys(R, F)
    
    for cand_k in cand_keys:
        if cand_k in keys:
            str_output += "Candidate key " + cand_k + " is a key\n"
        else:
            str_output += "Candidate key " + cand_k + " is NOT a key\n"
    
    return str_output

# function: return true if a schema is in 3NF
def is_3NF(R, F, keys):
    keys = find_keys(R, F)     # find keys 
    if not keys:    # if they still don't exist, return error
        print("Couldn't find any key. Exiting...")
        exit(1)
    else:
        print("Keys: " + str(keys))
    '''
    per ogni parte a sinistra
        - controllo se parte sinistra è una superchiave: ciclo ogni chiave e controllo se chiave in superchiave
            - se sì, continuo a ciclare nel principale, mi sta bene
            - altrimenti
                - per ogni parte a destra
                    - controllo se parte destra fa parte di una chiave: ciclo ogni chiave e vado a vedere se la parte a destra è nella singola chiave
                        - se sì, setto il booleano a true e continuo il ciclo principale
                        - se no, continuo questo ciclo qui

                        Se arrivo qua vuol dire che il ciclo è finito, quindi devo andare a controllare se il booleano dei primi è settato a false
                        se sì: esco
    '''
    # scorro ogni parte a sinistra
    for left_side, list_right_sides in F.items():
        # all'inizio entrambe le condizioni di 3NF sono a false
        superk = False
        prime = False

        if verbose : print("Checking " + left_side + "...")

        # controllo che a sinistra ci sia una superchiave 
        for key in keys:
            if key in left_side:
                superk = True
                break

        # visto che non è superchiave, devo controllare se gli attributi a destra sono primi
        if not superk:
            if verbose : print("\t" + left_side + " is not a superkey. Analyzing dependencies (primality test)...")
            
            # a noi ci basta trovare se una parte destra è contenuta in una delle chiavi, non ci interessa analizzarle tutte
            if not prime:   # IMPORTANTE: se non mettessimo questo controlleremmo anche tutte le altre chiavi, ma è inutile
                for right_side in list_right_sides:
                    if verbose : print("\n" + "\t" + "Analyzing: " + left_side + " -> " + right_side)
                    for key in keys:
                        if right_side in key:
                            if verbose: print("\t" + "\t" + right_side + " is in " + key)
                            prime = True
                            break       # IMPORTANTE: se non mettessimo questo controlleremmo anche tutte le altre chiavi, ma è inutile
                        else:
                            if verbose: print("\t" + "\t" + right_side + " is not in " + key)

                    # check se 3NF: se nessuno dei due è a True, lo schema non è per forza in 3NF, per definizione.
                    if not prime and not superk:
                        return False
        else:
            if verbose : print(left_side + " is a superkey. (contained in " + key + ")\n")

    return True

# fa il merge parte a destra e sinistra DA UN DIZIONARIO --> potenzialmente utile nell'algo di decomp
def merge_left_right(min_coverage):
    merged_list = list()

    # merge parti destre e sinistre
    for parte_destra in min_coverage:
        parte_sinistra = min_coverage[parte_destra]
        merged_list.append(set(parte_destra).union(set(parte_sinistra)))
    
    return merged_list

# fa il merge parte a destra e sinistra DA UNA STRINGA --> potenzialmente utile nell'algo di decomp
def merge_left_right_str(F):
    merged_list = list()

    # merge parti destre e sinistre
    final = re.split(', ', F)

    for df in final:
        nested_union = set()
        df = re.split(' -> ', df)
        for d in df: 
            nested_union = nested_union.union(set(d))
        merged_list.append(nested_union)
    return merged_list

# controlla se una copertura minimale F preserva le DF
# TODO: questa funzione è da completare
def check_decomposition(R, F, ro):
    min_coverage = minimal_coverage(R, F)
    new_F = dict_to_str(min_coverage)
    print("\nMinimal coverage: " + new_F)
    
    # trasformiamo ogni sottos di ro in un insieme
    final = []
    if len(ro) == 1:
        almost_final = re.split(', ', ro[0])
        for el in almost_final:
            final.append(set(el))
    else:
        for r in ro:    # altrimenti l'utente ha messo una lista
            final.append(set(r))

    # testo a turno le df: se qualcosa non è uno dei due sottoschemi ritorno False
    # scorro il dizionario e mi ricavo un set che poi andro a testare 
    for parte_sinistra in min_coverage:
        parte_destra = min_coverage[parte_sinistra]
        set_parte_sinistra = set(parte_sinistra)
        for p in parte_destra:
            set_parte_destra = set_parte_sinistra.union(set(p))

            contained = False
            for sottoschema in ro:
                # adesso controllo se contenuto in sottos di ro
                if set_parte_destra.issubset(sottoschema):
                    contained = True
                    if verbose : print(str(set_parte_destra) + " is contained in ro")
                    break

            if not contained:
                cl = []
                print(str(set_parte_destra) + " needs to be checked!")
                
                # TODO: algoritmo di decomposizione
          
    return True

# stampa la tabella della decomp
def print_table(tabella, decomp, str_R):
    # intestazione
    s = "  "
    for i in range(len(str_R)):
        s += "\t" + str_R[i]
    print(s)

    for i in range(0, len(str_R)):
        # costruisco la riga
        s = ""
        for el in tabella[i]:
            s += "\t" + str(el)
        
        print(decomp[i] + ":" + s)
        i += 1


##############
# ORG. FISICA
##############
# dato R, F e ro, controlla se c'è un join senza perdita
# TODO: questa funzione è da completare
def lossless_check(R, F, ro):
    str_R = ''.join(sorted(list(R)))

    # inizializzo la tabella a zero
    w, h = len(R), len(ro)
    tabella = [[None] * w for i in range(h)]
    decomp = list(sorted(list(ro)))
    
    # creo la famosa tabella
    for i in range(0,w):
        for j in range(0,h):
            # riempio le righe
            if str_R[j] in set(decomp[i]):
                tabella[i][j] = "A" + str(j+1)
            else:
                tabella[i][j] = "B" + str(i+1) + str(j+1)

            j += 1
        i += 1

    print_table(tabella, decomp, str_R)
    return True

def physical_isam(phyjsondir, v):
    kind = ""
    NR = 0
    R = 0
    K = 0
    CB = 0
    P = 0
    filled_space = 0
    floor_min = False
    kind, NR, R, K, CB, P, filled_space, floor_min = json_load(phyjsondir, v, True)
    
    if kind.upper() != "ISAM":
        print("The file you are trying does not contain an ISAM test.\nMake sure the file is correct and try again!")
        exit(1)

    MP = int(CB / R)
    MI = int(CB / (K + P))
    print("\nNR: \t" + str(NR) +
    "\nR: \t" + str(R) +
    "\nK: \t" + str(K) + 
    "\nCB: \t" + str(CB) +
    "\nP: \t" + str(P) +
    "\nFilled space: \t" + str(filled_space) +
    "\nFloor to min: \t" + str(floor_min))

    if verbose: print("MP: " + str(MP))
    if verbose: print("MI: " + str(MI))

    if filled_space == 100:
        BFP = math.ceil(NR/MP)
        BFI = math.ceil(BFP/MI)
    else:
        if not floor_min:
            CB_reduced = math.ceil((CB * filled_space /100))
            MP_reduced = math.ceil(CB_reduced / R)
            MI_reduced = math.ceil(CB_reduced/(K+P))
        else:
            CB_reduced = int((CB * filled_space /100))
            MP_reduced = int(CB_reduced / R)
            MI_reduced = int(CB_reduced/(K+P))
        
        MI = int(CB / (K + P))                

        BFP = math.ceil(NR/MP_reduced)
        BFI = math.ceil(BFP/MI_reduced)
    
    if not floor_min:
        MAX = math.ceil(math.log(BFI, 2)) + 1
    else:
        MAX = int(math.log(BFI, 2)) + 1
    
        if verbose: print("MP_reduced: " + str(MP_reduced))
        if verbose: print("MI_reduced: " + str(MI_reduced))
        if verbose: print("Filled space: " + str(filled_space) + "% --> CB" + str(filled_space) + ": " + str(CB_reduced))
    
    print("3A\tBFP \t# of blocks for main file:\t" + str(BFP))
    print("3B\tBFI \t# of blocks for index:\t\t" + str(BFI))
    print("3C\tMAX \tBinary search cost:\t\t" + str(MAX))

def physical_hash(phyjsondir, v):
    kind = "HASH"
    NR = 0
    R = 0
    K = 0
    CB = 0
    P = 0
    filled_space = 0
    floor_min = False
    kind, NR, R, K, CB, P, filled_space, floor_min, B, no_access = json_load(phyjsondir, v, True)
    
    if kind.upper() != "HASH":
        print("The file you are trying does not contain an ISAM test.\nMake sure the file is correct and try again!")
        exit(1)

    print("\nNR: \t" + str(NR) +
    "\nR: \t" + str(R) +
    "\nK: \t" + str(K) + 
    "\nCB: \t" + str(CB) +
    "\nP: \t" + str(P) +
    "\nFilled space: \t" + str(filled_space) +
    "\nFloor to min: \t" + str(floor_min) + 
    "\nB: \t" + str(B))

    RB = math.ceil(NR / B) 
    RBF = int((CB-P) / R)
    
    BB = math.ceil(RB / RBF)
    PB = int(CB / P)
    BBD = math.ceil(B / PB)

    MR = int(BB / 2)
    
    # blocksize bucket directory + total bucket size
    BLKTOT = BBD + (B * BB)

    if verbose: print("Number of record per bucket: (RB)\t\t" + str(RB))
    if verbose: print("Full records into a bucket: (RBF)\t\t" + str(RBF))
    if verbose: print("Block per bucket: (BB)\t\t" + str(BB))
    if verbose: print("Pointers per bucket: (PB)\t\t" + str(PB))
    if verbose: print("Block for bucket directory: (BBD)\t\t" + str(BBD))
    
    print("\n")
    print("3A\tBLKTOT \tSize of Bucket Directory and hash file: " + str(BLKTOT))
    print("3B\tMR \tSearch cost:\t\t" + str(MR))
    
    if no_access != 0:
        NBC = math.ceil(NR / (4 * (no_access * 2)))
        print("3C\tNBC \tNo bucket to have a search cost less equal than " + str(no_access) +": " + str(NBC))
    else:
        print("3C\tNBC \tNumber of accesses not provided.")

def main():
    global verbose
    global very_verbose
    global jsondir
    global phyjsondir
    global isam_opt
    global hash_opt

    # dichiarazione variabili 
    R = set()
    F = dict()
    X = list()
    keys = list()
    ro = list()
    all_files = False

    shortargs = "hvxn:k:t:c:m:p:l:"
    longargs = [
        'help',
        'verbose',
        'extended-verbose',
        'third-nf=',
        'third-fn=',
        '3FN=',
        '3NF=',
        'find-keys=',
        'test-keys=',
        'closure=',
        'minimal=',
        'json=',
        'preserve=',
        'lossless=',
        "physical-isam=",
        "physical-hash=",
        "physical-btree=",
        "phyjson=",
        "extended-help"
    ]

    # se vuoto
    if len(sys.argv[1:]) == 0:
        usage()

    # se tutti i file
    if "all" in sys.argv[1:]:
        all_files = True

    try:
        opts, args = getopt.getopt(sys.argv[1:], shortargs, longargs)
    except getopt.GetoptError as err:
        # print help
        print("Invalid option or argument given")
        usage()
        sys.exit(2)
        
    (opts, argv) = getopt.getopt(sys.argv[1:], shortargs, longargs)
    for (k, v) in opts:
        if k in ("-h", "--help"):
            usage()
        elif k in ("--extended-help"):
            usage(True)
        elif k in ("-v", "--verbose"):
            verbose = True
        elif k in ("-x", "--extended-verbose"):
            verbose = True
            very_verbose = True
        elif k in ("--json"):
            jsondir = v
        elif k in ("--phyjson"):
            phyjsondir = v
        elif k in ("-n", "--third-nf", "--third-fn", "--3FN", "--3NF"):
            print("CHECK 3NF")
            print("=========")
            if not all_files:
                R, F, keys, ro, X = json_load(jsondir, v)
                print_header(R, F, keys, ro, X)
                print("YES, this scheme is in 3NF" if is_3NF(R, F, keys) else "NO, this scheme IS NOT in 3NF")
            else:
                print("Checking if every scheme is in 3NF.")
                for f in os.listdir(jsondir):
                    if f.endswith(".json"):
                        os.system('cls' if os.name == 'nt' else 'clear')
                        print("Analyzing " + f)
                        R, F, keys, ro, X = json_load(jsondir, f)
                        print_header(R, F, keys, ro, X)
                        print("YES, this scheme is in 3NF" if is_3NF(R, F, keys) else "NO, this scheme IS NOT in 3NF")
                        print("=============================") 
                        input("Press Enter to continue...")
                        print("\n") 
        elif k in ("-k", "--find-keys"):
            print("FINDING KEYS")
            print("============")
            if not all_files:
                R, F, keys, ro, X = json_load(jsondir, v)
                print_header(R, F, keys, ro, X)
                print("Keys: " + str(find_keys(R, F)))
            else:
                print("Finding all keys all files.")
                for f in os.listdir(jsondir):
                    if f.endswith(".json"):
                        os.system('cls' if os.name == 'nt' else 'clear')
                        print("Finding keys in " + f)
                        R, F, keys, ro, X = json_load(jsondir, f)
                        print_header(R, F, keys, ro, X)
                        print("Keys: " + str(find_keys(R, F)))
                        print("=============================") 
                        input("Press Enter to continue...")
                        print("\n")
        elif k in ("-t", "--test-keys"):
            print("TESTING KEYS")
            print("============")
            if not all_files:
                R, F, keys, ro, X = json_load(jsondir, v)
                print_header(R, F, keys, ro, X)

                print_err(keys, "keys")
        
                str_out = are_keys(R, F, keys)
                if str_out == "": str_out = "No candidate key is key."
                print(str_out)

                print("=============================") 
                print("\n") 
            else:
                print("Figuring out if candidate keys are keys.")
                for f in os.listdir(jsondir):
                    if f.endswith(".json"):
                        os.system('cls' if os.name == 'nt' else 'clear')
                                    
                        print("TESTING KEYS")
                        print("============")
                    
                        print("Finding keys in " + f)
                        
                        R, F, keys, ro, X = json_load(jsondir, f)
                        print_header(R, F, keys, ro, X)

                        print_err(keys, "keys")        
                        
                        str_out = are_keys(R, F, keys)
                        if str_out == "": str_out = "No candidate key is key."
                        print(str_out)
                        print("=============================") 
                        input("Press Enter to continue...")
                        print("\n") 
        elif k in ("-c", "--closure"):
            print("CLOSURES")
            print("========")
            if not all_files:
                R, F, keys, ro, X = json_load(jsondir, v)
                print_header(R, F, keys, ro, X)
                print_err(X, "X")
                print("Closure: " + str(closures_calc(R, F, X)))
            else:
                print("Figuring out closures.")
                for f in os.listdir(jsondir):
                    if f.endswith(".json"):
                        os.system('cls' if os.name == 'nt' else 'clear')
                                                
                        print("CLOSURES")
                        print("========")
                    
                        print("Finding closures in " + f)
                        
                        R, F, keys, ro, X = json_load(jsondir, f)
                        print_header(R, F, keys, ro, X)

                        print_err(keys, "keys")        
                        
                        print("Closure: " + str(closures_calc(R, F, X)))
                        print("=============================") 
                        input("Press Enter to continue...")

        elif k in ("-m", "--minimal"):
            print("MINIMAL COVERAGE")
            print("================")
            if not all_files:
                R, F, keys, ro, X = json_load(jsondir, v)
                print_header(R, F, keys, ro, X)
                print("Minimal coverage: " + dict_to_str(minimal_coverage(R, F)))
            else:
                for f in os.listdir(jsondir):
                    if f.endswith(".json"):
                        os.system('cls' if os.name == 'nt' else 'clear')
                                                
                        print("MINIMAL COVERAGE")
                        print("================")
                        
                        R, F, keys, ro, X = json_load(jsondir, f)
                        print_header(R, F, keys, ro, X)
                               
                        print("Minimal coverage: " + dict_to_str(minimal_coverage(R, F)))
                        print("=============================") 
                        input("Press Enter to continue...")
        elif k in ("-p", "--preserve"):
            print("CHECK DECOMP")
            print("============")
            if not all_files:
                R, F, keys, ro, X = json_load(jsondir, v)
                print_header(R, F, keys, ro, X)
                # @TODO :commento cio: quando sara' implementato stampera' qualcosa di giusto, ripristinare questa riga
                # print("YES, ro preserve FDs" if check_decomposition(R, F, ro) else "NO, decompositions DO NOT preserve FDs")
            else:
                print("Checking if every scheme is in 3NF.")
                for f in os.listdir(jsondir):
                    if f.endswith(".json"):
                        os.system('cls' if os.name == 'nt' else 'clear')
                        print("Analyzing " + f)
                        R, F, keys, ro, X = json_load(jsondir, f)
                        print_header(R, F, keys, ro, X)
                        # @TODO :commento cio: quando sara' implementato stampera' qualcosa di giusto, ripristinare questa riga
                        # print("YES, the decompositions preserve FDs" if check_decomposition(R, F, ro) else "NO, decompositions DO NOT preserve FDs")
                        print("=============================") 
                        input("Press Enter to continue...")
                        print("\n") 
        elif k in ("-l", "--lossless"):
            print("CHECK LOSSLESS")
            print("==============")
            if not all_files:
                R, F, keys, ro, X = json_load(jsondir, v)
                print_header(R, F, keys, ro, X)
                
                # @TODO :commento cio: quando sara' implementato stampera' qualcosa di giusto, ripristinare questa riga
                #print("YES, lossless join is present" if lossless_check(R, F, ro) else "NO, lossless join IS NOT present")
            else:
                print("Checking if every scheme is in 3NF.")
                for f in os.listdir(jsondir):
                    if f.endswith(".json"):
                        os.system('cls' if os.name == 'nt' else 'clear')
                        print("Analyzing " + f)
                        R, F, keys, ro, X = json_load(jsondir, f)
                        print_header(R, F, keys, ro, X)
                        # @TODO :commento cio: quando sara' implementato stampera' qualcosa di giusto, ripristinare questa riga
                        #print("YES, lossless join is present" if lossless_check(R, F, ro) else "NO, lossless join IS NOT present")
                        print("=============================") 
                        input("Press Enter to continue...")
                        print("\n") 

        ######################
        # PHYSICAL EXERCISES #
        ######################
        elif k in ("--physical-isam"):
            isam_opt = True
            if all_files:        
                print("Calculating ISAM data for every JSON file")
                for f in os.listdir(phyjsondir):
                    if f.endswith(".json"):
                        os.system('cls' if os.name == 'nt' else 'clear')
                        physical_isam(phyjsondir, f)
                        input("Press Enter to continue...")
                        print("\n") 
            else:
                physical_isam(phyjsondir, v)
        elif k in ("--physical-hash"):
            hash_opt = True
            if all_files:        
                print("Calculating hash data for every JSON file")
                for f in os.listdir(phyjsondir):
                    if f.endswith(".json"):
                        os.system('cls' if os.name == 'nt' else 'clear')
                        physical_hash(phyjsondir, f)
                        input("Press Enter to continue...")
                        print("\n") 
            else:
                physical_hash(phyjsondir, v)            
        elif k in ("--physical-btree"):
            btree_opt = True
            print("coming soon!")
            pass

        else:
            assert False, "unhandled option"
    
if __name__ == "__main__":
    main()