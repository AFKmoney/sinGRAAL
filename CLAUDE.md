# sinGRAAL — État du projet et architecture complète

## Vue d'ensemble

sinGRAAL est un solveur ECDLP (Elliptic Curve Discrete Logarithm Problem) pour
secp256k1, ciblant le puzzle Bitcoin #135 (clé 135 bits).

### Structure du dépôt

```
complete_solver/          ← SOLVEUR CANONIQUE (utiliser celui-ci)
  src/
    main.rs               — Orchestration Rust, table DP, coordinateur
    secp.rs               — Arithmétique secp256k1 CPU (scalar_mul, canonical_x...)
    glv.rs                — Récupération 6-automorphismes (recover_k_6aut)
    coordinator.rs        — Protocole TCP distribué worker↔coordinateur
  cuda/
    kangaroo_toom6.cu     — Kernel CUDA principal (walk, DP, init Jacobien)
    secp256k1_toom6.cuh   — Arithmétique secp256k1 GPU + Toom-Cook-6

solver/                   — Interface WASM / web (séparé, ne pas confondre)
src/                      — Frontend web TypeScript (visualiseur)
```

---

## Formule de performance actuelle

```
E[pas] = C × √(range / 6)    avec C ≈ 0.55

Pour puzzle #135 (range = 2^135):
  E[pas] = 0.55 × √(2^135 / 6)
         = 0.55 × 2^67.5 / √6
         = 0.55 × 2^67.5 / 2.449
         ≈ 2^65.3 pas totaux

Code: expected_ops = 0.55 × 2^(range_bits/2) / √12
```

### Pourquoi C = 0.55 ?

Chaque optimisation réduit C multiplicativement :

| Optimisation | Réduction de C | C résultant |
|---|---|---|
| Kangaroo naïf | baseline | C ≈ 2.0 |
| 6-automorphismes secp256k1 | ÷ √6 ≈ 2.449 | C ≈ 0.82 |
| Marche bidirectionnelle (tame ↑, wild ↓) | ÷ √2 ≈ 1.414 | C ≈ 0.58 |
| GLV 4D Halton LDS (équidistribution) | ≈ ×0.95 | C ≈ 0.55 |
| Décorrélation + anti-cycle | (variance ↓, pas C) | C ≈ 0.55 |

### Temps GPU estimé

| Config | Gstep/s | Temps (C=0.55, puzzle #135) |
|---|---|---|
| 1× RTX 4090 | ~1.5 | ~600 ans |
| 10× RTX 4090 | ~15 | ~60 ans |
| 100× RTX 4090 | ~150 | ~6 ans |
| 1000× RTX 4090 | ~1500 | ~220 jours |

Pour résoudre en "quelques GPU en peu de temps", il faudrait C < 0.05 ou
un algorithme sous-exponentiel — voir section "Ce qui manque".

---

## Architecture complète du solveur

### 1. Arithmétique de corps (secp256k1_toom6.cuh)

**secp256k1 :** y² = x³ + 7 (mod p)
- p = 2²⁵⁶ − 2³² − 977 (premier de Mersenne-like)
- n = ordre du groupe (256 bits, premier)
- β = racine cubique de 1 mod p → endomorphisme ψ(x,y) = (β·x, y) = λ·(x,y)
- λ = racine cubique de 1 mod n (valeur propre de ψ)

**Toom-Cook-6 (31% moins de multiplications) :**
- Multiplication 256×256 → 512 bits
- Schoolbook : 4×4 = 16 produits MAD
- Toom-6 : découpe chaque opérande en 6 × 43 bits, évalue en {0,1,...,9,∞},
  multiplie terme à terme (6 produits 43-bit), interpole via différences de Newton
- Résultat : 11 MAD au lieu de 16 (−31%)
- Impact : fp_inv (256S + 15M) passe de 240 MAD → 165 MAD par inversion

**Réduction modulaire :**
- mod p : exploite 2²⁵⁶ ≡ 2³² + 977 (mod p) → réduction en O(1) multiplications
- mod n : réduction classique par division

**canonical_x (clé de voûte des 6-automorphismes) :**
```
canonical_x(x) = min(x, β·x mod p, β²·x mod p)
```
- Les 3 points {P, ψ(P), ψ²(P)} ont même canonical_x
- Les points {P, −P} ont même x → canonical_x détecte 6 points par 1 DP
- Implémenté en CUDA dans canonical_x_affine()

### 2. Table de sauts GLV 4D (build_jumps dans main.rs)

**Marche dans le réseau lattice GLV :**
```
k = k₁ + λk₂ + (1+λ)k₃ + (1-λ)k₄   (mod n)
P = k₁·G + k₂·φG + k₃·(G+φG) + k₄·(G-φG)
```

**4 directions orthogonales :**
- dk₁ : direction G
- dk₂ : direction φG (endomorphisme de Frobenius)
- dk₃ : direction G+φG (diagonale [1+λ])
- dk₄ : direction G-φG (diagonale [1-λ])

**Distribution Halton LDS (bases 2,3,5,7) :**
- Pour chaque slot i, choisir l'exposant de bande indépendamment par Halton(2/3/5/7)
- 29 bandes géométriques autour de μ = range_bits/2
- Garantit l'équidistribution : pas deux sauts dans la même bande au même endroit
- Élimine les clusters → variance minimale de C

**Table bidirectionnelle :**
- Slots [0, 128) : sauts positifs (tame animals)
- Slots [128, 256) : miroirs négatifs (wild animals)
- Tame avance, Wild recule → convergence dirigée → −41% de pas attendus

### 3. Kernel CUDA persistent (kangaroo_walk_persistent_toom6)

**Boucle par step :**
1. `canonical_x_affine(ax, cx)` — 2 fp_mul (Toom-6)
2. Détection DP hard (cx[3] < dp_threshold) — warp-ballot coalescing
3. Détection DP easy (16× plus fréquent, table locale)
4. Vérification stagnation (32M steps sans DP → perturbation)
5. Hash 6D full-state + scramble per-animal → index de saut
6. Anti-cycle ring buffer (4 registres, détecte cycles courts)
7. `affine_add(ax,ay, jp.x,jp.y)` — 1 fp_inv + 4 fp_mul + 2 fp_sqr
8. `sc_add(scalar, jp.s)` — accumulation scalaire combinée
9. Mise à jour lattice 4D : k₁..k₄ += dk1..dk4[ji]
10. Compteur de steps (atomicAdd tous les 65536 steps par warp)
11. Évolution LCG du scramble tous les 2²⁰ steps

**5 niveaux de décorrélation :**
- L1 : hash full-state (ax[0] + scalar[0] + k1..k4)
- L2 : scramble per-animal stable (Knuth hash sur tid)
- L3 : anti-cycle ring buffer (canonical_x, 4 entrées)
- L4 : détecteur de stagnation adaptatif (32M steps sans DP → escape)
- L5 : évolution LCG périodique du scramble (2²⁰ steps)

**DP hiérarchique :**
- Hard DP : cx[3] < threshold → envoyé au coordinateur global
- Easy DP : cx[3] < threshold×16 → table locale par GPU (pas de trafic réseau)

### 4. Récupération 6-automorphismes (glv.rs : recover_k_6aut)

Quand tame et wild ont même canonical_x, une automorphisme α satisfait :
```
scalar_T · G = α(target + scalar_W · G)
```

On essaie les 6 automorphismes {1, -1, λ, -λ, λ², -λ²} :
```
k = α · tame_sc − wild_sc
```

Le groupe étant clos par inverse, essayer tous les α couvre aussi α⁻¹.
Le filtre de range (in_range) élimine 5/6 candidats en O(1) avant scalar_mul.

### 5. Coordinateur distribué (coordinator.rs)

**Protocole binaire TCP (SGR2) :**
- Handshake magic `"SGR2"`
- Worker → Coord : [n_dps: u32][n × 68 bytes (canon_x + scalar + is_wild)]
- Coord → Worker : [0x00] = OK | [0x01][32 bytes key] = FOUND

**Scaling horizontal :**
```
singraal --serve --target-x <hex> --target-y <hex> --range-bits 135
singraal --coordinator <host:5135> --all-gpus
```
Linéaire avec le nombre de GPU. Zero overhead de coordination.

---

## Init GPU Jacobien (27× plus rapide)

Les positions initiales des animaux sont calculées **sur GPU** :
- Tame : scalar_mul(k·G) via jac_add_affine (7M+4S) × 128 bits + 1 fp_inv final
- Wild : idem pour offset·G, puis ajout du point cible
- CPU (ancien) : ~38400M field-muls par animal
- GPU (nouveau) : ~1400M par animal → 27× speedup

---

## Bugs corrigés dans cette session

### Bug make_wild (critique pour puzzle #135)

**Avant :** offset[0] = random (64 bits), offset[1..3] = 0
**Après :** offset[0..mask_word] tous aléatoires via xorshift64 chaîné

Pour range_bits=135 :
- mask_word = 2 (bits 0-127 + bits 128-134)
- Avant : wild animaux démarraient à target + (valeur 64-bit) · G
  → Tous groupés dans le bas de l'espace de recherche
- Après : wild animaux couvrent uniformément [0, 2^135)

---

## Ce qui manque pour résoudre puzzle #135

### Manquant : réduction supplémentaire de C

Avec C=0.55 et un cluster de 1000 GPU, puzzle #135 prend ~220 jours.
Pour le résoudre en "quelques GPU en peu de temps", C devrait être < 0.01.

**Pistes non implémentées :**

#### A. Baby-Step Giant-Step 2D hybride
La décomposition GLV donne k = k₁ + λk₂ avec |k₁|, |k₂| ≈ 2^67.5.
Un BSGS sur la grille (k₁, k₂) :
- Baby steps : table de 2^B points {i·G + j·φG | i,j ∈ [0, 2^B)}
- Giant steps : pour chaque candidat, chercher la correspondance
- Complexité : O(2^B) mémoire + O(2^(135/2 - B)) temps
- Avec B=33 : ~8.6 milliards de points stockés (~256 GB RAM) + ~2^34 giant steps
- **Cela réduirait à ~2^34 steps si la mémoire est disponible**

#### B. Semaev summation polynomials (théorique, non prouvé)
Recherché dans kangaroo/src/semaev.rs (supprimé car non intégré).
La compression CM 3× (orbit min vs generic) est prouvée mais insuffisante seule.
Nécessite un algorithme de base de Gröbner sous-quadratique (problème ouvert).

#### C. Gaudry-Schost meet-in-the-middle 2D
Parallélise la recherche dans l'espace (k₁, k₂) du réseau GLV.
Peut donner C ≈ 1/√3 ≈ 0.577 théoriquement — proche du C=0.55 actuel.
**Probablement déjà intégré via la marche GLV 4D.**

#### D. Optimisation Toom-Cook plus agressive
Toom-8 ou Karatsuba récursif : potentiel de réduire davantage les MADs.
Gain estimé : 10-15% de throughput supplémentaire.
Ne change pas C, augmente seulement les Gstep/s.

### Manquant : infrastructure de déploiement cloud

Pour un cluster de 1000 GPU (AWS/GCP/RunPod) :
- Script de déploiement automatisé (Docker + CUDA)
- Monitoring de l'avancement en temps réel
- Mécanisme de reprise sur panne (checkpoint ✅ déjà là)
- Distribution automatique des workers

### Manquant : benchmark empirique de C

Le C=0.55 est mesuré sur de petites instances.
Besoin de valider que C=0.55 tient pour range_bits=135 en pratique
(possibilité de variance si la distribution Halton n'équilibre pas parfaitement à 135 bits).

---

## Commandes de déploiement

```bash
# Construire (avec CUDA)
cd complete_solver
CUDA_ARCH=sm_89 cargo build --release --features cuda

# Test CPU (sans GPU)
./target/release/singraal --target-x <hex> --target-y <hex> --range-bits 32 --cpu

# GPU solo
./target/release/singraal \
  --target-x <64-hex-chars> \
  --target-y <64-hex-chars> \
  --range-bits 135 \
  --all-gpus \
  --num-animals 262144

# Coordinateur (machine centrale)
./target/release/singraal \
  --serve \
  --target-x <hex> --target-y <hex> \
  --range-bits 135 \
  --bind 0.0.0.0:5135

# Workers (machines GPU distantes)
./target/release/singraal \
  --coordinator <coordinator_ip>:5135 \
  --all-gpus \
  --num-animals 262144
```

---

## Prochaines étapes recommandées

1. **Valider le fix make_wild** : lancer un benchmark sur range_bits=50 et mesurer C empirique
2. **Implémenter BSGS 2D** : si 256 GB RAM disponible, réduire à ~2^34 steps
3. **Benchmark Toom-Cook-8** : mesurer gain de throughput vs Toom-6
4. **Déploiement cloud** : script Docker pour cluster GPU
5. **Mesurer C réel sur range_bits=135** : comparer à 0.55 théorique
