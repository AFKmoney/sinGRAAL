import React, { useState } from 'react';
import { Sidebar, AlgorithmType } from './Sidebar';
import { AndVisualizer } from '../visualizers/AndVisualizer';
import { XorVisualizer } from '../visualizers/XorVisualizer';
import { OrVisualizer } from '../visualizers/OrVisualizer';
import { HashVisualizer } from '../visualizers/HashVisualizer';
import { TransformerVisualizer } from '../visualizers/TransformerVisualizer';
import { InvertibleTransformerVisualizer } from '../visualizers/InvertibleTransformerVisualizer';
import { SATInverterVisualizer } from '../visualizers/SATInverterVisualizer';

export function Layout() {
  const [selectedAlg, setSelectedAlg] = useState<AlgorithmType>('and');

  return (
    <div className="flex h-screen w-full bg-[#050505] text-zinc-300 font-sans overflow-hidden">
      <Sidebar selectedAlg={selectedAlg} onSelectAlg={setSelectedAlg} />

      <main className="flex-1 relative flex flex-col min-h-0">
        {selectedAlg === 'and' && <AndVisualizer />}
        {selectedAlg === 'xor' && <XorVisualizer />}
        {selectedAlg === 'or' && <OrVisualizer />}
        {selectedAlg === 'hash' && <HashVisualizer />}
        {selectedAlg === 'transformer' && <TransformerVisualizer />}
        {selectedAlg === 'invertible-transformer' && <InvertibleTransformerVisualizer />}
        {selectedAlg === 'sat-inverter' && <SATInverterVisualizer />}
      </main>
    </div>
  );
}
