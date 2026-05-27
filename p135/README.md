# sinGRAAL — Solveur Bitcoin Puzzle #135

Solveur ECDLP secp256k1 production-ready pour le Bitcoin Puzzle #135 (clé 135 bits).
Pipeline algébrique complet : GLV + Semaev S₃ + LLL Kill-Switch certifié + GSDD + 6-automorphismes.

---

## Architecture du pipeline

```
k·G = P    (k ∈ [2^134, 2^135))
```

### Étape 1 — Décomposition GLV 4D
```
k = k₁ + λk₂ + (1+λ)k₃ + (1-λ)k₄  (mod n)
```
- λ = endomorphisme secp256k1 (racine cubique de 1 mod n)
- Réduit l'espace de recherche 2D : |k₁|, |k₂| ≈ 2^(range/2)
- 4 directions orthogonales couvrent le réseau GLV complet

### Étape 2 — Polynôme de Semaev S₃
Pour k = A + δ + λ(B + ε), le polynôme S₃ satisfait :
```
S₃(x(A·G + δ·G), x(B·φG + ε·φG), x_P) = 0  (mod p)
```
Développé en série bivariée en (δ, ε) → matrice de Macaulay **28×28** (niveau m=3, Jochemsz-May).

### Étape 3 — LLL Kill-Switch certifié (innovation clé)
Critère de rejet certifié :
```
λ₁(L)² ≥ p⁶/28  →  tile TUÉE  (aucune racine dans [A, A+X) × [B, B+X))
```
- Early-abort LLL : arrêt dès que `min‖b*_k‖² ≥ borne` (O(µs) par tile)
- Kill-rate cible : >99.99% des tiles éliminées sans brute-force
- Tiles survivantes : brute-force O(X²) avec X = 2^block_bits

### Étape 4 — 6-automorphismes secp256k1
Les 6 automorphismes `{±1, ±λ, ±λ²}` agissent sur P :
```
φ(x, y) = (β·x, y)  →  φ(P) = λ·P
```
- Chaque tile couvre 6 cibles `{P, -P, λP, -λP, λ²P, -λ²P}` simultanément
- ×6 couverture de l'espace k par tile → ×6 moins de tiles à visiter

### Étape 5 — GSDD (Galois Symmetry + Nested Field Decomposition)
- Exposant de Frobenius : `d = p mod (n-1)` (d ≈ 2^129 pour secp256k1)
- Décomposition CRT sur les petits facteurs de (n-1) : {2, 3, 149, 631, ...}
- Cantor-Zassenhaus (Tonelli-Shanks) pour les racines de polynômes mod n

### Étape 6 — AnchorTable L2
- Pré-calcul de tous les points d'ancrage `ia·step·G` et `ib·step·φG`
- Tient en cache L1/L2 CPU → zéro scalar_mul pendant la recherche principale

### Étape 7 — Rayon parallel
- Boucle externe sur `ia` distribuée sur tous les cœurs CPU
- Work-stealing automatique (Rayon) → scaling linéaire avec N cœurs

---

## Prérequis

| Composant | Version minimale |
|---|---|
| Docker | 20.10+ |
| Docker Compose | v2.0+ |
| RAM | 16 GB (block_bits=20), 64 GB (block_bits=24) |
| CPU | 8+ cœurs recommandés |
| CUDA (optionnel) | 12.0+ pour le filtre GPU |

---

## Déploiement rapide

### 1. Obtenir les coordonnées du puzzle #135

Les coordonnées du point cible Bitcoin Puzzle #135 sont publiques :
```bash
# Remplacer par les vraies valeurs du puzzle
export TARGET_X=<coordonnée_x_en_hex_64_chars>
export TARGET_Y=<coordonnée_y_en_hex_64_chars>
```

### 2. Construire l'image Docker

```bash
# Depuis la racine du dépôt
cd /chemin/vers/sinGRAAL

docker build \
  -f p135/Dockerfile \
  -t singraal-p135:latest \
  .
```

### 3. Lancer le solveur

```bash
docker run --rm \
  -e TARGET_X=$TARGET_X \
  -e TARGET_Y=$TARGET_Y \
  -e RANGE_BITS=135 \
  -e BLOCK_BITS=20 \
  -e THREADS=0 \
  -v $(pwd)/results:/data \
  singraal-p135:latest
```

### 4. Vérifier avec selftest (avant de lancer sur #135)

```bash
# Test automatique sur clé 40 bits (instantané)
docker run --rm -e TARGET_X=selftest singraal-p135:latest

# Selftest GSDD complet
docker run --rm -e TARGET_X=gsdd-selftest singraal-p135:latest
```

---

## Déploiement cloud GPU

### RunPod

```bash
# Lancer N pods GPU (ex: RTX 4090)
runpodctl create pod \
  --name "singraal-p135" \
  --imageName "votre-registry/singraal-p135:latest" \
  --gpuType "NVIDIA GeForce RTX 4090" \
  --containerDiskSize 20 \
  --env "TARGET_X=$TARGET_X" \
  --env "TARGET_Y=$TARGET_Y" \
  --env "RANGE_BITS=135" \
  --env "BLOCK_BITS=20"
```

### vast.ai

```bash
# Chercher instances GPU disponibles
vastai search offers 'gpu_name=RTX_4090 num_gpus=1 reliability>0.98'

# Louer et déployer
vastai create instance <OFFER_ID> \
  --image "votre-registry/singraal-p135:latest" \
  --env "-e TARGET_X=$TARGET_X -e TARGET_Y=$TARGET_Y"
```

### Lambda Labs

```bash
# Dans le terminal de l'instance Lambda
git clone https://github.com/AFKmoney/sinGRAAL
cd sinGRAAL
docker build -f p135/Dockerfile -t singraal-p135 .
docker run -d \
  -e TARGET_X=$TARGET_X \
  -e TARGET_Y=$TARGET_Y \
  -v /home/ubuntu/results:/data \
  singraal-p135:latest
```

### Multi-instances (docker-compose)

```bash
# Lancer 4 instances en parallèle
cd sinGRAAL
TARGET_X=$TARGET_X TARGET_Y=$TARGET_Y \
  docker compose -f p135/docker-compose.yml up --scale solver=4 -d

# Surveiller les logs
docker compose -f p135/docker-compose.yml logs -f
```

---

## Variables d'environnement

| Variable | Défaut | Description |
|---|---|---|
| `TARGET_X` | *(requis)* | Coordonnée x du point cible (hex 64 chars) |
| `TARGET_Y` | *(requis)* | Coordonnée y du point cible (hex 64 chars) |
| `RANGE_BITS` | `135` | Taille de l'espace de recherche (k < 2^RANGE_BITS) |
| `BLOCK_BITS` | `20` | Taille des tiles (X = 2^BLOCK_BITS). RAM : ~16 GB à 20, ~256 GB à 24 |
| `THREADS` | `0` | Nombre de threads CPU (0 = tous les cœurs disponibles) |

---

## Résultats attendus

### Sortie standard

```
╔═══════════════════════════════════════════════════════════╗
║  FULL STACK — toutes innovations actives                  ║
╠═══════════════════════════════════════════════════════════╣
║  #1 GLV-4D  #2 Semaev-S₃  #3 Frobenius-CRT              ║
║  #4 LLL-m=3-kill  #5 6-aut  #6 L2-anchor  #7 Rayon      ║
╠═══════════════════════════════════════════════════════════╣
║  range=135  block=20  half=70  tiles≈2^...
╚═══════════════════════════════════════════════════════════╝
[full-stack] anchors prêts en 0.XX s
[full-stack] 0.01%  kill=99.9X%  surv=N  t=XX.Xs
...
[full-stack] ✓ SURVIVANT a=... b=... norm²≈2^XX  t=XXXs
FOUND k = 0x<clé_hex>
```

### Fichier solution

En cas de succès, la clé est sauvegardée dans `/data/solution.txt` :
```
# sinGRAAL — Bitcoin Puzzle #135 Solution
# Date: ...
TARGET_X=...
TARGET_Y=...
RANGE_BITS=135
k=<clé_hex>
```

---

## Paramétrage BLOCK_BITS

| BLOCK_BITS | Taille tile X | RAM (anchors) | Kill-rate LLL | Temps par tile survivante |
|---|---|---|---|---|
| 16 | 65 536 | ~1 MB | très élevé | ~1 ms |
| 20 | 1 048 576 | ~16 MB | élevé | ~1 s |
| 22 | 4 194 304 | ~64 MB | moyen | ~16 s |
| 24 | 16 777 216 | ~256 MB | variable | ~256 s |

**Recommandation** : commencer avec `BLOCK_BITS=20` pour mesurer le kill-rate réel,
puis ajuster selon les résultats.

---

## 10 Innovations mathématiques

| # | Innovation | Description |
|---|---|---|
| 1 | Index Calculus artificiel | Base de facteurs sur la courbe pour décomposer les points |
| 2 | Semaev S₃ | Polynôme trivarié S₃(x₁,x₂,x₃)=0 iff Σ points = O |
| 3 | Gröbner F4/F5 | Résolution du système polynomial par bases de Gröbner parallélisées |
| 4 | Weil Descent virtuel | Descente sur corps premiers via extensions Fp^k |
| 5 | GLV 4D endomorphisme | k = k₁+λk₂+(1+λ)k₃+(1-λ)k₄, 4 directions orthogonales |
| 6 | LLL m=3 kill-switch | Matrice Macaulay 28×28, rejet certifié λ₁²≥p⁶/28 |
| 7 | Block Lanczos / Wiedemann | Algèbre linéaire creuse sur tenseur secp256k1 |
| 8 | Contraintes CRT sémantiques | Décomposition k mod qᵢ via BSGS sur E[qᵢ] |
| 9 | Réduction par isogénies | Transfert vers courbes isogènes plus faibles |
| 10 | GSDD Frobenius | kᵈ≡k (mod n), d=p mod(n-1), contrainte algébrique supplémentaire |

---

## Licence

Propriété de Philippe-Antoine Robert. Usage académique et de recherche.
