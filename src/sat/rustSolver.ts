// Bridge between TypeScript encoders and the Rust WASM CDCL solver.
// Provides the same API surface as the TypeScript solveCDCL() but backed by Rust.

import { WasmSolver } from '../wasm/solver';
import type { SolveStats, SolveResult } from './sha256sat';

export type { SolveStats as RustSolveStats };

// Pack clause array into flat Int32Array with 0 separators
function packClauses(clauses: number[][]): Int32Array {
  let size = 0;
  for (const c of clauses) size += c.length + 1;
  const buf = new Int32Array(size);
  let pos = 0;
  for (const c of clauses) {
    for (const l of c) buf[pos++] = l;
    buf[pos++] = 0;
  }
  return buf;
}

// Main solver function — same signature as solveCDCL from sha256sat.ts
export async function solveRust(
  allClauses: number[][],
  numVars: number,
  onProgress?: (stats: SolveStats, assigned: number) => void,
  cancelRef?: { cancelled: boolean },
): Promise<SolveResult> {
  const solver = new WasmSolver(numVars);
  const packed = packClauses(allClauses);
  solver.add_clauses_packed(packed);

  const STEPS_PER_FRAME = 2000; // CDCL conflict budget per animation frame

  return new Promise((resolve) => {
    let rafHandle = 0;

    function tick() {
      if (cancelRef?.cancelled) {
        solver.free();
        resolve({ status: 'cancelled', assignment: new Map(), stats: parseStats('{}') });
        return;
      }

      const json = solver.step(STEPS_PER_FRAME);
      const data = JSON.parse(json) as {
        status: string;
        decisions: number;
        propagations: number;
        conflicts: number;
        learned: number;
        backjumps: number;
        maxBackjump: number;
        assigned: number;
        numVars: number;
      };

      const stats: SolveStats = {
        decisions: data.decisions,
        propagations: data.propagations,
        conflicts: data.conflicts,
        learnedClauses: data.learned,
        backjumps: data.backjumps,
        maxBackjump: data.maxBackjump,
      };

      if (onProgress) onProgress(stats, data.assigned);

      if (data.status === 'sat') {
        const rawAsgn = solver.get_assignment();
        const assignment = new Map<number, boolean>();
        for (let i = 0; i < rawAsgn.length; i++) {
          if (rawAsgn[i] !== 0) assignment.set(i + 1, rawAsgn[i] === 1);
        }
        solver.free();
        resolve({ status: 'sat', assignment, stats });
      } else if (data.status === 'unsat') {
        solver.free();
        resolve({ status: 'unsat', assignment: new Map(), stats });
      } else {
        rafHandle = requestAnimationFrame(tick);
      }
    }

    rafHandle = requestAnimationFrame(tick);
  });
}

function parseStats(json: string): SolveStats {
  try {
    const d = JSON.parse(json);
    return {
      decisions: d.decisions ?? 0,
      propagations: d.propagations ?? 0,
      conflicts: d.conflicts ?? 0,
      learnedClauses: d.learned ?? 0,
      backjumps: d.backjumps ?? 0,
      maxBackjump: d.maxBackjump ?? 0,
    };
  } catch { return { decisions:0, propagations:0, conflicts:0, learnedClauses:0, backjumps:0, maxBackjump:0 }; }
}
