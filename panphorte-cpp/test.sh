#!/bin/bash

# Configurazione
CPP_SOURCE="main1.cpp"
CPP_EXECUTABLE="panphorte"
PYTHON_SCRIPT="main_walks.py"
INPUT_GFA="../chr9B_walks_fix.gfa"
# INPUT_GFA="../example3_walks.gfa"

# Colori per output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=== Automazione Test Python e C++ ==="
echo ""

# 1. Compila il C++ solo se modificato o se l'eseguibile non esiste
if [ ! -f "$CPP_EXECUTABLE" ] || [ "$CPP_SOURCE" -nt "$CPP_EXECUTABLE" ]; then
    echo -e "${YELLOW}Compilazione di $CPP_SOURCE...${NC}"
    g++ -std=c++17 -Wall -g -o "$CPP_EXECUTABLE" "$CPP_SOURCE"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Compilazione completata${NC}"
    else
        echo -e "${RED}✗ Errore di compilazione${NC}"
        exit 1
    fi
else
    echo -e "${GREEN}✓ Eseguibile C++ già aggiornato${NC}"
fi

echo ""

# 2. Esegui il programma Python
# echo -e "${YELLOW}Esecuzione di $PYTHON_SCRIPT...${NC}"
# echo "--- Output Python ---"
# /usr/bin/time -f "Time (Python): %e %E, MAX Memory: %M KB" python3 ../"$PYTHON_SCRIPT" -i "$INPUT_GFA"
# PYTHON_EXIT=$?

# if [ $PYTHON_EXIT -eq 0 ]; then
#     echo -e "${GREEN}✓ Python terminato con successo${NC}"
# else
#     echo -e "${RED}✗ Python terminato con errore (codice: $PYTHON_EXIT)${NC}"
# fi

# echo ""

# 3. Esegui il programma C++
echo -e "${YELLOW}Esecuzione di $CPP_EXECUTABLE...${NC}"
echo "--- Output C++ ---"
/usr/bin/time -f "Time (C++): %e %E, MAX Memory: %M KB" ./"$CPP_EXECUTABLE" -i "$INPUT_GFA"
CPP_EXIT=$?

if [ $CPP_EXIT -eq 0 ]; then
    echo -e "${GREEN}✓ C++ terminato con successo${NC}"
else
    echo -e "${RED}✗ C++ terminato con errore (codice: $CPP_EXIT)${NC}"
fi

echo ""
echo "=== Test completato ==="