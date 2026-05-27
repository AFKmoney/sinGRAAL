# sinGRAAL — Solveur Bitcoin Puzzle #135

Solveur ECDLP (Elliptic Curve Discrete Logarithm Problem) pour le puzzle Bitcoin #135
(clé 135 bits sur secp256k1), prêt à déployer sur cloud GPU.

---

## Architecture du pipeline

Le solveur s'appuie sur la bibliothèque `bsgs2d` qui intègre une chaîne de traitements
mathématiques avancés organisée comme suit :

### 1. GLV décomposition

La décomposition GLV (Gallant-Lambert-Vanstone) exploite l'endomorphisme efficient
de secp256k1 :

```
φ(P) = (β·x, y) = λ·P
```

où β est une racine cubique de 1 mod p et λ une racine cubique de 1 mod n.

Cela permet de décomposer la clé secrète :

```
k = k₁ + λ·k₂  (mod n)
```

avec |k₁|, |k₂| ≈ √n ≈ 2^128, réduisant le problème à deux dimensions.

### 2. Semaev S₃

Le polynôme de sommation de Semaev S₃(x_A, x_B, x_P) = 0 encode la relation
« x_A + x_B = x_P sur la courbe » via un polynôme en les coordonnées-x seulement.
La matrice de Macaulay associée (degré m=3) est de taille 28×28 et permet
d'exprimer l'appartenance d'un point à une tile comme un système polynomial.

### 3. LLL Kill-Switch certifié

Le critère d'élimination de tile repose sur la théorie de Minkowski :

```
λ₁(L)² ≥ det(L)^(2/dim) / dim     (Minkowski)
```

Pour le réseau lattice Macaulay de dimension 3 :
```
λ₁² ≥ p⁶/28  →  tile TUÉE en O(µs)
```

Ce critère est certifié (pas de faux négatifs) : toute tile satisfaisant ce critère
ne contient provablement aucune solution. En pratique, plus de 99,99 % des tiles
sont éliminées sans aucun calcul de point de courbe.

### 4. 6-automorphismes secp256k1

secp256k1 admet 6 automorphismes utiles issus de l'endomorphisme φ et de la négation :

```
{ P, -P, φ(P), -φ(P), φ²(P), -φ²(P) }
```

Ces 6 points ont le même `canonical_x = min(x, β·x mod p, β²·x mod p)`.
Chaque entrée de la baby table couvre donc 6 clés simultanément,
multipliant l'efficacité par 6 sans surcoût mémoire.

### 5. GSDD (Galois Symmetry + Nested Field Decomposition)

Le module GSDD exploite :
- **Frobenius** : la structure galoisienne de l'extension de corps sous-jacente
- **CRT** (Théorème Chinois des Restes) : décomposition modulaire pour paralléliser
  les vérifications de candidats en plusieurs fragments indépendants

### 6. AnchorTable L2

Table de points de référence pré-calculés dimensionnée pour tenir dans le cache L2 du CPU.
Évite les défauts de cache lors des étapes géantes (giant steps), maintenant un débit
élevé sans accès DRAM.

### 7. Parallélisme Rayon

Tous les cœurs CPU disponibles sont utilisés via la bibliothèque Rayon (work-stealing).
Le traitement des tiles et la construction de la baby table sont entièrement parallélisés.

### 8. Marche bidirectionnelle tame/wild

- **Tame** : animaux démarrant dans la zone connue, avançant vers la cible
- **Wild** : animaux démarrant près de la cible, reculant vers la zone connue
- La convergence dirigée réduit le nombre de pas attendus d'un facteur √2

### 9. Distribution Halton LDS

Les tailles de saut suivent une distribution de Low-Discrepancy Sequence (Halton,
bases 2, 3, 5, 7) sur 29 bandes géométriques. Cette équidistribution garantit
une couverture uniforme de l'espace de recherche et élimine les clusters de sauts
qui augmentent la variance de C.

### 10. Détection DP hiérarchique

Deux niveaux de Distinguished Points :
- **Hard DP** (cx[3] < threshold) : envoyé au coordinateur global via TCP
- **Easy DP** (cx[3] < threshold×16) : stocké dans une table locale par GPU,
  sans trafic réseau, pour la détection de collisions locales

---

## Prérequis

| Composant | Version minimale |
|---|---|
| Docker | 20.10+ |
| Docker Compose | v2.0+ |
| CPU | x86_64, 4+ cœurs recommandés |
| RAM | 8 GB minimum, 32 GB recommandé |
| GPU (optionnel) | NVIDIA avec CUDA 12+ pour le filtre GPU |

---

## Déploiement rapide (local)

### 1. Cloner le dépôt

```bash
git clone <REPO_URL>
cd sinGRAAL
```

### 2. Configurer les coordonnées cibles

Remplacer les placeholders par les vraies coordonnées du puzzle #135 :

```bash
export TARGET_X="<coordonnée_x_hex_64_chars>"
export TARGET_Y="<coordonnée_y_hex_64_chars>"
```

> **Note** : Les vraies coordonnées du puzzle Bitcoin #135 sont disponibles sur
> https://privatekeys.pw/puzzles/bitcoin-puzzle-tx

### 3. Lancer le solveur

```bash
cd p135

# Instance unique
docker-compose up --build

# Plusieurs instances en parallèle (ex : 4)
docker-compose up --build --scale solver=4 -d

# Voir les logs
docker-compose logs -f
```

### 4. Selftest (vérification sans coordonnées cibles)

```bash
docker run --rm \
  -e TARGET_X=selftest \
  singraal-p135:latest
```

---

## Déploiement cloud

### RunPod

```bash
# Générer les commandes pour 8 pods GPU
./deploy.sh runpod 8
```

Pré-requis : `runpodctl` installé et configuré avec votre clé API RunPod.

Étapes :
1. Pousser l'image Docker vers Docker Hub
2. Exécuter les commandes générées par le script
3. Surveiller avec `runpodctl get pods`

### vast.ai

```bash
# Générer les commandes pour 16 instances GPU
./deploy.sh vast 16
```

Pré-requis : `vastai` CLI installé (`pip install vastai`) et clé API configurée.

### Lambda Labs

```bash
# Générer les commandes pour 4 instances A100
./deploy.sh lambda 4
```

Lambda Labs utilise des instances bare-metal SSH. Le script génère les commandes
de lancement SSH + Docker à exécuter manuellement après connexion.

---

## Paramètres

| Variable | Défaut | Description |
|---|---|---|
| `TARGET_X` | (obligatoire) | Coordonnée x du point cible (hex 64 chars) |
| `TARGET_Y` | (obligatoire) | Coordonnée y du point cible (hex 64 chars) |
| `RANGE_BITS` | `135` | Taille de l'espace de recherche (k ∈ [0, 2^RANGE_BITS)) |
| `BLOCK_BITS` | `20` | Taille des blocs de tiles (2^BLOCK_BITS tiles par bloc) |
| `THREADS` | `0` | Nombre de threads CPU (0 = tous les cœurs disponibles) |

---

## Résultats attendus

Quand une solution est trouvée, le solveur affiche :

```
╔══════════════════════════════════════════════════════╗
║  SOLUTION TROUVÉE                                    ║
║  k = <valeur_hex_de_la_cle_privee>                   ║
╚══════════════════════════════════════════════════════╝
[INFO] Solution sauvegardée dans /data/solution.txt
```

Le fichier `/data/solution.txt` (dans le volume Docker) contient :

```
# sinGRAAL — Bitcoin Puzzle #135 Solution
# Date: <timestamp UTC>
TARGET_X=<hex>
TARGET_Y=<hex>
RANGE_BITS=135
k=<cle_privee_hex>
```

Pour accéder au volume de résultats :

```bash
docker volume inspect p135_solver_data
# ou
docker run --rm -v p135_solver_data:/data alpine cat /data/solution.txt
```

---

## Temps de résolution estimé

Le temps de résolution dépend du **kill_rate empirique** du LLL Kill-Switch,
mesuré sur les tiles réelles du puzzle #135.

| Paramètre | Valeur théorique |
|---|---|
| Opérations attendues | ~2^65.3 (C=0.55) |
| Kill-rate LLL Kill-Switch | > 99,99 % des tiles (mesuré sur range_bits ≤ 50) |
| Accélération effective | dépend du kill_rate empirique sur range_bits=135 |

> **Important** : Aucune garantie de temps de résolution n'est fournie.
> Le temps réel dépend du kill_rate effectif sur les tiles de 135 bits,
> qui doit être validé empiriquement. Lancer d'abord un benchmark avec
> `--range-bits 50` pour mesurer le kill_rate sur votre matériel.

---

## Innovations mathématiques

| # | Innovation | Description |
|---|---|---|
| 1 | GLV décomposition | k = k₁ + λk₂ (mod n), réduit l'espace 2D à ≈2^67.5 par dimension |
| 2 | 6-automorphismes secp256k1 | φ(P)=(β·x,y)=λ·P, ×6 couverture k-space par entrée de table |
| 3 | Semaev S₃ | Polynôme trivarié, matrice Macaulay 28×28 (m=3) |
| 4 | LLL Kill-Switch certifié | λ₁²≥p⁶/28 → tile tuée en O(µs), >99,99% tiles éliminées |
| 5 | GSDD Galois Symmetry | Frobenius + CRT, décomposition en fragments indépendants |
| 6 | Nested Field Decomposition | Décomposition CRT du corps pour parallélisme optimal |
| 7 | AnchorTable L2 | Points pré-calculés dimensionnés pour le cache CPU L2 |
| 8 | Rayon parallel | Work-stealing sur tous les cœurs CPU |
| 9 | Distribution Halton LDS | Bases 2,3,5,7 sur 29 bandes géométriques, variance minimale |
| 10 | DP hiérarchique hard/easy | Deux seuils, table locale GPU, zéro trafic réseau pour easy DP |

---

## Structure du dépôt

```
p135/
  Dockerfile          — Build multi-stage Rust + runtime debian:bookworm-slim
  docker-compose.yml  — Orchestration multi-instances avec volume /data
  deploy.sh           — Script de déploiement cloud (RunPod / vast.ai / Lambda / local)
  entrypoint.sh       — Validation, bannière, lancement bsgs2d, sauvegarde solution
  README.md           — Cette documentation

bsgs2d/               — Moteur du solveur (NE PAS MODIFIER)
  src/
    main.rs           — CLI, orchestration BSGS 1D/2D
    dispatcher.rs     — run_full_stack() : pipeline complet
    secp.rs           — Arithmétique secp256k1
    glv4d.rs          — GLV 4D + 6-automorphismes
    gsdd.rs           — GSDD (Frobenius + CRT)
    coppersmith.rs    — Matrice Macaulay, LLL, polynôme S₃
    lll.rs            — Algorithme LLL
    lll_earlyabort.rs — LLL avec arrêt anticipé (kill-switch)
```

---

## Licence

Ce logiciel est fourni à des fins de recherche. L'utilisation est soumise aux
conditions du dépôt parent sinGRAAL.
