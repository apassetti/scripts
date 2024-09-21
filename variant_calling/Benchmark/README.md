# Benchmark Pipeline

Questa pipeline è relativa alla creazione di read sintetiche utilizzando il programma DWGSIM. Le read sintetiche sono basate sul genoma di riferimento hg38 e specificamente sul cromosoma 7.

## Requisiti

- DWGSIM
- hg38 reference genome

## Istruzioni

1. Scaricare e installare DWGSIM.
2. Scaricare il genoma di riferimento hg38.
3. Eseguire DWGSIM specificando il cromosoma 7 come regione di interesse.

## Esempio di comando

```bash
dwgsim -1 100 -2 100 -r 0.001 -R 0.15 -X 0.3 -y 0.02 -z 100 -N 1000000 hg38_chr7.fasta chr7_reads
```

## Specifica dei comandi

- `-1 100`: Lunghezza delle read del primo file (forward) in bp.
- `-2 100`: Lunghezza delle read del secondo file (reverse) in bp.
- `-r 0.001`: Tasso di errore di base.
- `-R 0.15`: Tasso di errore di lettura.
- `-X 0.3`: Tasso di errore di inserimento.
- `-y 0.02`: Tasso di errore di delezione.
- `-z 100`: Seed per il generatore di numeri casuali.
- `-N 1000000`: Numero di read da generare.
- `hg38_chr7.fasta`: File del genoma di riferimento.
- `chr7_reads`: Prefisso per i file di output delle read generate.


## Output

Il comando genererà file di read sintetiche che possono essere utilizzati per ulteriori analisi di benchmark.

## Contatti

Per ulteriori informazioni, contattare [Andrea_Passetti] all'indirizzo [andrea.passetti91@gmail.com].

# Pipeline

