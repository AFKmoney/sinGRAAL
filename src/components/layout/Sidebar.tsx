import React from 'react';
import { ArrowLeftRight, Split, GitCommit, Network, Filter, Layers, ScanSearch, KeyRound, Cpu, Rabbit } from 'lucide-react';

export type AlgorithmType = 'and' | 'or' | 'xor' | 'hash' | 'transformer' | 'invertible-transformer' | 'sat-inverter' | 'ec-puzzle' | 'hybrid-solver' | 'kangaroo';

interface SidebarProps {
  selectedAlg: AlgorithmType;
  onSelectAlg: (alg: AlgorithmType) => void;
}

const navItems: { id: AlgorithmType; label: string; icon: React.ReactNode; desc: string }[] = [
  { id: 'and', label: 'AND Gate', icon: <Filter className="w-4 h-4" />, desc: 'Lossy Information' },
  { id: 'or', label: 'OR Gate', icon: <Split className="w-4 h-4" />, desc: 'Lossy Information' },
  { id: 'xor', label: 'XOR Gate', icon: <ArrowLeftRight className="w-4 h-4" />, desc: 'Perfectly Reversible' },
  { id: 'hash', label: 'SHA-256 Hash', icon: <GitCommit className="w-4 h-4" />, desc: 'Pseudo-Chaos Pipeline' },
  { id: 'transformer', label: 'Transformer (ReLU)', icon: <Network className="w-4 h-4" />, desc: 'Non-Linear Loss' },
  { id: 'invertible-transformer', label: 'Invertible Transformer', icon: <Layers className="w-4 h-4" />, desc: 'XOR-based Reversal Matrix' },
  { id: 'sat-inverter', label: 'SHA-256 Réduit', icon: <ScanSearch className="w-4 h-4" />, desc: 'Recherche de Préimage' },
  { id: 'ec-puzzle', label: 'EC Puzzle Solver', icon: <KeyRound className="w-4 h-4" />, desc: 'ECDLP via CNF/CDCL' },
  { id: 'hybrid-solver', label: 'Hybrid SAT+BSGS', icon: <Cpu className="w-4 h-4" />, desc: '135 → 2×2⁴⁰ Complexity' },
  { id: 'kangaroo', label: 'Kangaroo 6-Aut', icon: <Rabbit className="w-4 h-4" />, desc: '√6 speedup via ψ²' },
];

export function Sidebar({ selectedAlg, onSelectAlg }: SidebarProps) {
  return (
    <aside className="w-72 bg-zinc-900/20 border-r border-zinc-800 flex flex-col h-full shrink-0">
      <div className="p-6 border-b border-zinc-800">
        <h1 className="text-xl font-light tracking-widest text-white uppercase">
          Aetherial Inversion <span className="text-[10px] text-zinc-500 font-mono align-top ml-1 mt-1 inline-block">v4.0.2</span>
        </h1>
        <p className="text-[10px] text-zinc-500 mt-1 uppercase tracking-tighter">Visualizer Matrix</p>
      </div>

      <nav className="flex-1 overflow-y-auto p-4 space-y-2">
        {navItems.map((item) => {
          const isSelected = selectedAlg === item.id;
          return (
            <button
              key={item.id}
              onClick={() => onSelectAlg(item.id)}
              className={`w-full text-left p-4 rounded-lg transition-all duration-200 border ${
                isSelected
                  ? 'bg-zinc-900/50 border-cyan-500/30'
                  : 'bg-transparent border-transparent hover:bg-zinc-900/30 hover:border-zinc-700/50'
              }`}
            >
              <div className="flex items-center space-x-3 mb-1">
                <div className={`p-1.5 rounded text-[10px] uppercase font-bold tracking-widest ${isSelected ? 'text-cyan-400' : 'text-zinc-500'}`}>
                  {item.icon}
                </div>
                <span className={`text-xs uppercase font-bold tracking-widest ${isSelected ? 'text-zinc-100' : 'text-zinc-400'}`}>
                  {item.label}
                </span>
              </div>
              <p className="text-[9px] uppercase tracking-widest text-zinc-600 ml-11">{item.desc}</p>
            </button>
          );
        })}
      </nav>

      <div className="p-6 border-t border-zinc-800">
        <div className="bg-zinc-900/50 rounded flex flex-col p-4 border border-zinc-800">
          <div className="text-[10px] text-zinc-600 uppercase font-bold tracking-widest mb-1">Quantum State</div>
          <div className="flex items-center space-x-2">
            <div className="w-1.5 h-1.5 rounded-full bg-emerald-500 animate-pulse shadow-[0_0_8px_rgba(16,185,129,0.5)]"></div>
            <span className="text-[10px] text-emerald-400 font-mono uppercase tracking-widest">Coherent</span>
          </div>
        </div>
      </div>
    </aside>
  );
}
