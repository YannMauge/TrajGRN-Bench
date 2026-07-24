import numpy as np
from .gene import gene
from scipy.stats import ttest_rel, ttest_ind, ranksums
import sys
import csv
import networkx as nx
from scipy.stats import wasserstein_distance

class sergio (object):

    def __init__(self,number_genes, number_bins, number_sc, noise_params,\
    noise_type, decays, dynamics = False, sampling_state = 10, tol = 1e-3,\
    window_length = 100, dt = 0.01, optimize_sampling = False,\
    bifurcation_matrix = None, noise_params_splice = None, noise_type_splice = None,\
    splice_ratio = 4, dt_splice = 0.01, migration_rate = None,\
    max_memory_mb = 4096.0):
        """
        Noise is a gaussian white noise process with zero mean and finite variance.
        noise_params: The amplitude of noise in CLE. This can be a scalar to use
        for all genes or an array with the same size as number_genes.
        Tol: p-Value threshold above which convergence is reached
        window_length: length of non-overlapping window (# time-steps) that is used to realize convergence
        dt: time step used in  CLE
        noise_params and decays: Could be an array of length number_genes, or single value to use the same value for all genes
        number_sc: number of single cells for which expression is simulated
        sampling_state (>=1): single cells are sampled from sampling_state * number_sc steady-state steps
        optimize_sampling: useful for very large graphs. If set True, may help finding a more optimal sampling_state and so may ignore the input sampling_state
        noise_type: We consider three types of noise, 'sp': a single intrinsic noise is associated to production process, 'spd': a single intrinsic noise is associated to both
        production and decay processes, 'dpd': two independent intrinsic noises are associated to production and decay processes
        dynamics: whether simulate splicing or not
        bifurcation_matrix: is a numpy array (nBins_ * nBins) of <1 values; bifurcation_matrix[i,j] indicates whether cell type i differentiates to type j or not. Its value indicates the rate of transition. If dynamics == True, this matrix should be specified
        noise_params_splice: Same as "noise_params" but for splicing. if not specified, the same noise params as pre-mRNA is used
        noise_type_splice: Same as "noise_type" but for splicing. if not specified, the same noise type as pre-mRNA is used
        splice_ratio: it shows the relative amount of spliced mRNA to pre-mRNA (at steady-state) and therefore tunes the decay rate of spliced mRNA as a function of unspliced mRNA. Could be an array of length number_genes, or single value to use the same value for all genes
        dt_splice = time step for integrating splice SDE
        max_memory_mb: per-gene-bin concentration buffer limit in MiB (default 4096 = 4 GiB).
            Raises MemoryError before allocation if any buffer would exceed this.


        Note1: It's assumed that no two or more bins differentiate into the same new bin i.e. every bin has either 0 or 1 parent bin
        Note2: differentitation rates (e.g. type1 -> type2) specified in bifurcation_matrix specifies the percentage of cells of type2 that are at the vicinity of type1
        """

        self.max_memory_mb_ = max_memory_mb
        self.nGenes_ = number_genes
        self.nBins_ = number_bins
        self.nSC_ = number_sc
        self.sampling_state_ = sampling_state
        self.tol_ = tol
        self.winLen_ = window_length
        self.dt_ = np.float32(dt)
        self.optimize_sampling_ = optimize_sampling
        self.level2verts_ = {}
        self.gID_to_level_and_idx = {} # This dictionary gives the level and idx in self.level2verts_ of a given gene ID
        self.binDict = {} # This maps bin ID to list of gene objects in that bin; only used for dynamics simulations
        self.maxLevels_ = 0
        self.init_concs_ = np.zeros((number_genes, number_bins), dtype=np.float32)
        self.meanExpression = -1 * np.ones((number_genes, number_bins), dtype=np.float32)
        self.noiseType_ = noise_type
        self.dyn_ = dynamics
        self.nConvSteps = np.zeros(number_bins) # This holds the number of simulated steps till convergence
        if dynamics:
            self.bifurcationMat_ = np.array(bifurcation_matrix)
            self.binOrders_ = []
            self.binDict = {}
            for b in range(self.nBins_):
                self.binDict[b] = np.zeros(self.nGenes_,).tolist()
        ############
        # This graph stores for each vertex: parameters(interaction
        # parameters for non-master regulators and production rates for master
        # regulators), tragets, regulators and level
        ############
        self.graph_ = {}

        if np.isscalar(noise_params):
            self.noiseParamsVector_ = np.repeat(np.float32(noise_params), number_genes)
        elif np.shape(noise_params)[0] == number_genes:
            self.noiseParamsVector_ = np.asarray(noise_params, dtype=np.float32)
        else:
            print ("Error: expect one noise parameter per gene")


        if np.isscalar(decays) == 1:
            self.decayVector_ = np.repeat(np.float32(decays), number_genes)
        elif np.shape(decays)[0] == number_genes:
            self.decayVector_ = np.asarray(decays, dtype=np.float32)
        else:
            print ("Error: expect one decay parameter per gene")
            sys.exit()


        if self.dyn_:
            if (self.bifurcationMat_ == None).any():
                print ("Error: Bifurcation Matrix is missing")
                sys.exit()

            if noise_type_splice == None:
                self.noiseTypeSp_ = noise_type
            else:
                self.noiseTypeSp_ = noise_type_splice


            if dt_splice == None:
                self.dtSp_ = np.copy(self.dt_)
            else:
                self.dtSp_ = dt_splice


            if noise_params_splice == None:
                self.noiseParamsVectorSp_ = np.copy(self.noiseParamsVector_)
            elif np.isscalar(noise_params_splice):
                self.noiseParamsVectorSp_ = np.repeat(noise_params_splice, number_genes)
            elif np.shape(noise_params_splice)[0] == number_genes:
                self.noiseParamsVectorSp_ = noise_params_splice
            else:
                print ("Error: expect one splicing noise parameter per gene")
                sys.exit()

            if np.isscalar(splice_ratio):
                self.ratioSp_ = np.repeat(splice_ratio, number_genes)
            elif np.shape(splice_ratio)[0] == number_genes:
                self.ratioSp_ = splice_ratio
            else:
                print ("Error: expect one splicing ratio parameter per gene")
                sys.exit()





    def build_graph (self, input_file_taregts, input_file_regs, shared_coop_state = 0):
        """
        # 1- shared_coop_state: if >0 then all interactions are modeled with that
        # coop state, and coop_states in input_file_taregts are ignored. Otherwise,
        # coop states are read from input file. Reasonbale values ~ 1-3
        # 2- input_file_taregts: a csv file, one row per targets. Columns: Target Idx, #regulators,
        # regIdx1,...,regIdx(#regs), K1,...,K(#regs), coop_state1,...,
        # coop_state(#regs)
        # 3- input_file_regs: a csv file, one row per master regulators. Columns: Master regulator Idx,
        # production_rate1,...,productions_rate(#bins)
        # 4- input_file_taregts should not contain any line for master regulators
        # 5- For now, assume that nodes in graph are either master regulator or
        # target. In other words, there should not be any node with no incomming
        # or outgoing edge! OTHERWISE IT CAUSES ERROR IN CODE.
        # 6- The indexing of genes start from 0. Also, the indexing used in
        # input files should match the indexing (if applicable) used for initilizing
        # the object.
        """

        for i in range(self.nGenes_):
            self.graph_[i] = {}
            self.graph_[i]['targets'] = []


        allRegs = []
        allTargets = []

        with open(input_file_taregts,'r') as f:
            reader = csv.reader(f, delimiter=',')
            if (shared_coop_state <= 0):
                for row in reader:
                    nRegs = np.int(row[1])
                    ##################### Raise Error ##########################
                    if nRegs == 0:
                        print ("Error: a master regulator (#Regs = 0) appeared in input")
                        sys.exit()
                        ############################################################

                    currInteraction = []
                    currParents = []
                    for regId, K, C_state in zip(row[2: 2 + nRegs], row[2+nRegs : 2+2*nRegs], row[2+2*nRegs : 2+3*nRegs]):
                        currInteraction.append((np.int(regId), np.float(K), np.float(C_state), 0)) # last zero shows half-response, it is modified in another method
                        allRegs.append(np.int(regId))
                        currParents.append(np.int(regId))
                        self.graph_[np.int(regId)]['targets'].append(np.int(row[0]))

                    self.graph_[np.int(row[0])]['params'] = currInteraction
                    self.graph_[np.int(row[0])]['regs'] = currParents
                    self.graph_[np.int(row[0])]['level'] = -1 # will be modified later
                    allTargets.append(np.int(row[0]))

                    #if self.dyn_:
                    #    for b in range(self.nBins_):
                    #        binDict[b].append(gene(np.int(row[0]),'T', b))
            else:
                for indRow, row in enumerate(reader):
                    nRegs = np.int(np.float(row[1]))
                    ##################### Raise Error ##########################
                    if nRegs == 0:
                        print ("Error: a master regulator (#Regs = 0) appeared in input")
                        sys.exit()
                        ############################################################

                    currInteraction = []
                    currParents = []
                    for regId, K, in zip(row[2: 2 + nRegs], row[2+nRegs : 2+2*nRegs]):
                        currInteraction.append((np.int(np.float(regId)), np.float(K), shared_coop_state, 0)) # last zero shows half-response, it is modified in another method
                        allRegs.append(np.int(np.float(regId)))
                        currParents.append(np.int(np.float(regId)))
                        self.graph_[np.int(np.float(regId))]['targets'].append(np.int(np.float(row[0])))

                    self.graph_[np.int(np.float(row[0]))]['params'] = currInteraction
                    self.graph_[np.int(np.float(row[0]))]['regs'] = currParents
                    self.graph_[np.int(np.float(row[0]))]['level'] = -1 # will be modified later
                    allTargets.append(np.int(np.float(row[0])))

                    #if self.dyn_:
                    #    for b in range(self.nBins_):
                    #        binDict[b].append(gene(np.int(row[0]),'T', b))

        #self.master_regulators_idx_ = set(np.setdiff1d(allRegs, allTargets))

        with open(input_file_regs,'r') as f:
            masterRegs = []
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                if np.shape(row)[0] != self.nBins_ + 1:
                    print ("Error: Inconsistent number of bins")
                    sys.exit()

                masterRegs.append(int(float(row[0])))
                self.graph_[int(float(row[0]))]['rates'] = [np.float(i) for i in row[1:]]
                self.graph_[int(float(row[0]))]['regs'] = []
                self.graph_[int(float(row[0]))]['level'] = -1

                #if self.dyn_:
                #    for b in range(self.nBins_):
                #        binDict[b].append(gene(np.int(row[0]),'MR', b))

        self.master_regulators_idx_ = set(masterRegs)


        if (len(self.master_regulators_idx_) + np.shape(allTargets)[0] != self.nGenes_):
            print ("Error: Inconsistent number of genes")
            sys.exit()

        self.find_levels_(self.graph_) # make sure that this modifies the graph

        if self.dyn_:
            self.find_bin_order_(self.bifurcationMat_)

    def build_graph_from_rows(self, targets_rows, regs_rows, shared_coop_state=0):
        """
        In-memory variant of ``build_graph`` that accepts parsed rows
        (lists of lists) instead of CSV file paths.  Same semantics otherwise.
        """
        for i in range(self.nGenes_):
            self.graph_[i] = {}
            self.graph_[i]['targets'] = []

        allRegs = []
        allTargets = []

        if shared_coop_state <= 0:
            for row in targets_rows:
                nRegs = int(row[1])
                if nRegs == 0:
                    print("Error: a master regulator (#Regs = 0) appeared in input")
                    sys.exit()

                currInteraction = []
                currParents = []
                for regId, K, C_state in zip(
                    row[2: 2 + nRegs],
                    row[2 + nRegs: 2 + 2 * nRegs],
                    row[2 + 2 * nRegs: 2 + 3 * nRegs],
                ):
                    currInteraction.append((int(regId), float(K), float(C_state), 0))
                    allRegs.append(int(regId))
                    currParents.append(int(regId))
                    self.graph_[int(regId)]['targets'].append(int(row[0]))

                self.graph_[int(row[0])]['params'] = currInteraction
                self.graph_[int(row[0])]['regs'] = currParents
                self.graph_[int(row[0])]['level'] = -1
                allTargets.append(int(row[0]))
        else:
            for row in targets_rows:
                nRegs = int(float(row[1]))
                if nRegs == 0:
                    print("Error: a master regulator (#Regs = 0) appeared in input")
                    sys.exit()

                currInteraction = []
                currParents = []
                for regId, K in zip(
                    row[2: 2 + nRegs],
                    row[2 + nRegs: 2 + 2 * nRegs],
                ):
                    currInteraction.append(
                        (int(float(regId)), float(K), shared_coop_state, 0)
                    )
                    allRegs.append(int(float(regId)))
                    currParents.append(int(float(regId)))
                    self.graph_[int(float(regId))]['targets'].append(
                        int(float(row[0]))
                    )

                self.graph_[int(float(row[0]))]['params'] = currInteraction
                self.graph_[int(float(row[0]))]['regs'] = currParents
                self.graph_[int(float(row[0]))]['level'] = -1
                allTargets.append(int(float(row[0])))

        masterRegs = []
        for row in regs_rows:
            if len(row) != self.nBins_ + 1:
                print("Error: Inconsistent number of bins")
                sys.exit()

            masterRegs.append(int(float(row[0])))
            self.graph_[int(float(row[0]))]['rates'] = [float(i) for i in row[1:]]
            self.graph_[int(float(row[0]))]['regs'] = []
            self.graph_[int(float(row[0]))]['level'] = -1

        self.master_regulators_idx_ = set(masterRegs)

        if len(self.master_regulators_idx_) + len(allTargets) != self.nGenes_:
            print("Error: Inconsistent number of genes")
            sys.exit()

        self.find_levels_(self.graph_)

        if self.dyn_:
            self.find_bin_order_(self.bifurcationMat_)

    def find_levels_ (self, graph):
        """
        # This is a helper function that takes a graph and assigns layer to all
        # verticies. It uses longest path layering algorithm from
        # Hierarchical Graph Drawing by Healy and Nikolovself. A bottom-up
        # approach is implemented to optimize simulator run-time. Layer zero is
        # the last layer for which expression are simulated
        # U: verticies with an assigned layer
        # Z: vertizies assigned to a layer below the current layer
        # V: set of all verticies (genes)

        This also sets a dictionary that maps a level to a matrix (in form of python list)
        of all genes in that level versus all bins
        """

        U = set()
        Z = set()
        V = set(graph.keys())

        currLayer = 0
        self.level2verts_[currLayer] = []
        idx = 0

        # Pre-compute non-self targets for layering (self-loops for
        # auto-regulation must not block topological ordering).
        _non_self_targets = {}
        for v in V:
            _non_self_targets[v] = set(graph[v]['targets']) - {v}

        # Safety: detect genuine cycles (not just self-loops)
        _max_layers = len(V) + 1  # one more than possible layers
        while U != V:
            if currLayer > _max_layers:
                remaining = V - U
                raise RuntimeError(
                    f"find_levels_: cycle detected in GRN. "
                    f"Genes {sorted(remaining)} could not be layered. "
                    f"Check for mutual-regulation cycles (e.g. A→B and B→A)."
                )

            currVerts = set(
                filter(lambda v: _non_self_targets[v].issubset(Z), V - U)
            )

            # If stuck (mutual-regulation cycle), promote ALL remaining
            # vertices to the current layer so they are simulated together.
            # Within-cycle regulator reads will use same-step values, which
            # is a valid approximation for the Euler-Maruyama scheme.
            if not currVerts and U != V:
                remaining = list(V - U)
                remaining.sort(
                    key=lambda v: (
                        sum(1 for r in set(graph[v].get('regs', [])) if r != v and r in Z),
                        v,
                    ),
                    reverse=True,
                )
                currVerts = set(remaining)

            for v in currVerts:
                graph[v]['level'] = currLayer
                U.add(v)
                if {v}.issubset(self.master_regulators_idx_):
                    allBinList = [gene(v,'MR', i) for i in range(self.nBins_)]
                    self.level2verts_[currLayer].append(allBinList)
                    self.gID_to_level_and_idx[v] = (currLayer, idx)
                    idx += 1
                else:
                    allBinList = [gene(v,'T', i) for i in range(self.nBins_)]
                    self.level2verts_[currLayer].append(allBinList)
                    self.gID_to_level_and_idx[v] = (currLayer, idx)
                    idx += 1

            currLayer += 1
            Z = Z.union(U)
            self.level2verts_[currLayer] = []
            idx = 0

        self.level2verts_.pop(currLayer)
        self.maxLevels_ = currLayer - 1

        if not self.dyn_:
            self.set_scIndices_()

    def set_scIndices_ (self, safety_steps = 0):
        """
        # First updates sampling_state_ if optimize_sampling_ is set True: to optimize run time,
        run for less than 30,000 steps in first level
        # Set the single cell indices that are sampled from steady-state steps
        # Note that sampling should be performed from the end of Concentration list
        # Note that this method should run after building graph(and layering) and should
        be run just once!
        """

        if self.optimize_sampling_:
            state = np.true_divide(30000 - safety_steps * self.maxLevels_, self.nSC_)
            if state < self.sampling_state_:
                self.sampling_state_ = state

        self.scIndices_ = np.random.randint(low = - self.sampling_state_ * self.nSC_, high = 0, size = self.nSC_)

    def calculate_required_steps_(self, level, safety_steps = 0):
        """
        # Calculates the number of required simulation steps after convergence at each level.
        # safety_steps: estimated number of steps required to reach convergence (same), although it is not neede!
        """
        #TODO: remove this safety step

        return self.sampling_state_ * self.nSC_ + level * safety_steps

    def calculate_half_response_(self, level):
        """
        Calculates the half response for all interactions between previous layer
        and current layer.

        When a regulator has not been simulated yet (e.g. due to
        mutual-regulation cycles), a default half-response of 1.0 is used as
        a fallback so the simulation can proceed.
        """

        currGenes = self.level2verts_[level]

        for g in currGenes: # g is list of all bins for a single gene
            c = 0
            if g[0].Type == 'T':
                for interTuple in self.graph_[g[0].ID]['params']:
                    regIdx = interTuple[0]
                    meanArr = self.meanExpression[regIdx]

                    if set(meanArr) == set([-1]):
                        # Regulator not yet simulated (mutual-regulation
                        # cycle).  Use a default half-response so the
                        # simulation can proceed.
                        self.graph_[g[0].ID]['params'][c] = (
                            self.graph_[g[0].ID]['params'][c][0],
                            self.graph_[g[0].ID]['params'][c][1],
                            self.graph_[g[0].ID]['params'][c][2],
                            np.float32(1.0),
                        )
                    else:
                        self.graph_[g[0].ID]['params'][c] = (
                            self.graph_[g[0].ID]['params'][c][0],
                            self.graph_[g[0].ID]['params'][c][1],
                            self.graph_[g[0].ID]['params'][c][2],
                            np.float32(np.mean(meanArr)),
                        )
                    c += 1
            #Else: g is a master regulator and does not need half response

    def _estimate_step_memory_mb(self, n_steps, n_genes_level, n_bins=None):
        """Estimate memory (MiB) for concentration buffers at one level."""
        if n_bins is None:
            n_bins = self.nBins_
        return n_genes_level * n_bins * n_steps * 4 / (1024 * 1024)

    def hill_(self, reg_conc, half_response, coop_state, repressive=False, hr_cs=None):
        """
        Hill function for activation/repression.

        Accepts scalar or array ``reg_conc``.  Returns the same shape as
        ``reg_conc`` (float32).

        Parameters
        ----------
        hr_cs : float, optional
            Pre-computed ``half_response ** coop_state``.  When provided,
            avoids recomputing the power on every call.
        """
        reg_conc = np.asarray(reg_conc, dtype=np.float32)
        cs = np.float32(coop_state)

        if hr_cs is None:
            hr_cs = np.float32(np.power(np.float32(half_response), cs))

        # For reg_conc == 0: activation → 0, repression → 1
        zero_mask = (reg_conc == 0)
        pow_reg = np.power(reg_conc, cs)
        denom = np.float32(hr_cs) + pow_reg
        # Avoid division by zero
        denom = np.where(denom == 0, np.float32(1.0), denom)
        result = np.true_divide(pow_reg, denom)
        result = np.where(zero_mask, np.float32(0.0), result)

        if repressive:
            return np.float32(1.0) - result
        return result

    def init_gene_bin_conc_ (self, level):
        """
        Initilizes the concentration of all genes in the input level

        Note: calculate_half_response_ should be run before this method
        """

        currGenes = self.level2verts_[level]
        for g in currGenes:
            if g[0].Type == 'MR':
                allBinRates = self.graph_[g[0].ID]['rates']

                for bIdx, rate in enumerate(allBinRates):
                    g[bIdx].append_Conc(np.float32(np.true_divide(rate, self.decayVector_[g[0].ID])))

            else:
                params = self.graph_[g[0].ID]['params']

                for bIdx in range(self.nBins_):
                    rate = np.float32(0.0)
                    for interTuple in params:
                        meanExp = self.meanExpression[interTuple[0], bIdx]
                        rate += np.float32(np.abs(interTuple[1])) * self.hill_(meanExp, interTuple[3], interTuple[2], interTuple[1] < 0)

                    g[bIdx].append_Conc(np.float32(np.true_divide(rate, self.decayVector_[g[0].ID])))

    def calculate_prod_rate_(self, bin_list, level):
        """
        calculates production rates for the input list of gene objects in different bins but all associated to a single gene ID

        Returns a float32 numpy array of shape (n_bins,).
        """
        type = bin_list[0].Type

        if (type == 'MR'):
            rates = self.graph_[bin_list[0].ID]['rates']
            return np.array([rates[gb.binID] for gb in bin_list], dtype=np.float32)

        else:
            params = self.graph_[bin_list[0].ID]['params']
            Ks = np.array([np.abs(t[1]) for t in params], dtype=np.float32)
            regIndices = [t[0] for t in params]
            binIndices = [gb.binID for gb in bin_list]
            currStep = bin_list[0].simulatedSteps_
            lastLayerGenes = np.copy(self.level2verts_[level + 1])
            n_regs = len(regIndices)
            n_bins = len(binIndices)

            # --- Pre-gather all regulator concentrations ---
            reg_concs = np.empty((n_regs, n_bins), dtype=np.float32)
            for tupleIdx, rIdx in enumerate(regIndices):
                regGeneLevel = self.gID_to_level_and_idx[rIdx][0]
                regGeneIdx = self.gID_to_level_and_idx[rIdx][1]
                regGene_allBins = self.level2verts_[regGeneLevel][regGeneIdx]
                for colIdx, bIdx in enumerate(binIndices):
                    if regGene_allBins[bIdx]._use_buffer:
                        reg_concs[tupleIdx, colIdx] = regGene_allBins[bIdx].Conc[regGene_allBins[bIdx]._conc_ptr - 1]
                    else:
                        reg_concs[tupleIdx, colIdx] = regGene_allBins[bIdx].Conc[currStep]

            # --- Vectorized hill for all regulator×bin pairs ---
            hillMatrix = np.empty((n_regs, n_bins), dtype=np.float32)
            for tupleIdx in range(n_regs):
                hillMatrix[tupleIdx, :] = self.hill_(
                    reg_concs[tupleIdx, :],
                    params[tupleIdx][3],
                    params[tupleIdx][2],
                    params[tupleIdx][1] < 0,
                )

            return np.matmul(Ks, hillMatrix)


    @staticmethod
    def _compute_noise_(rng, rates, decay, noise_type, noise_param, sqrt_dt):
        """
        Compute stochastic noise term for the CLE.

        Parameters
        ----------
        rng : numpy.random.Generator
        rates : ndarray (..., n_bins) – production (or any additive) rates
        decay : ndarray (..., n_bins) – decay rates
        noise_type : str – 'sp', 'spd', or 'dpd'
        noise_param : float32 – noise amplitude
        sqrt_dt : float32 – sqrt(dt)

        Returns
        -------
        noise : ndarray, same shape as rates (already multiplied by sqrt_dt)
        """
        # Guard against tiny negative values from floating-point imprecision
        # before taking sqrt (mathematically rates/decay are >= 0).
        safe_rates = np.maximum(rates, np.float32(0.0))
        safe_decay = np.maximum(decay, np.float32(0.0))

        if noise_type == 'sp':
            dw = rng.normal(size=rates.shape).astype(np.float32)
            noise = noise_param * np.sqrt(safe_rates, dtype=np.float32) * dw
        elif noise_type == 'spd':
            dw = rng.normal(size=rates.shape).astype(np.float32)
            amplitude = noise_param * (
                np.sqrt(safe_rates, dtype=np.float32) + np.sqrt(safe_decay, dtype=np.float32)
            )
            noise = amplitude * dw
        else:  # 'dpd'
            dw_p = rng.normal(size=rates.shape).astype(np.float32)
            dw_d = rng.normal(size=rates.shape).astype(np.float32)
            noise = (
                noise_param * np.sqrt(safe_rates, dtype=np.float32) * dw_p
                + noise_param * np.sqrt(safe_decay, dtype=np.float32) * dw_d
            )
        return sqrt_dt * noise

    def CLE_simulator_(self, level):

        self.calculate_half_response_(level)
        nReqSteps = self.calculate_required_steps_(level)

        # ---- Pre-allocate concentration buffers BEFORE init (so initial conc
        #      is written into the pre-allocated buffer, not a Python list) ----
        sim_genes = self.level2verts_[level]
        n_genes_level = len(sim_genes)

        # --- Memory safety check ---
        est_mb = self._estimate_step_memory_mb(nReqSteps, n_genes_level)
        per_obj_mb = nReqSteps * 4 / (1024 * 1024)
        if est_mb > self.max_memory_mb_:
            raise MemoryError(
                f"Level {level}: estimated {est_mb:.1f} MiB for concentration "
                f"buffers exceeds max_memory_mb={self.max_memory_mb_:.0f} MiB. "
                f"({n_genes_level} genes × {self.nBins_} bins × {nReqSteps} steps). "
                f"Reduce sampling_state, n_cells, or increase max_memory_mb."
            )
        for g_list in sim_genes:
            for gb in g_list:
                gb.init_conc_buffer(nReqSteps, max_buffer_mb=per_obj_mb * 2)

        # ---- Now initialize concentrations into the pre-allocated buffer ----
        self.init_gene_bin_conc_(level)

        # ---- Pre-compute positive sc indices ----
        pos_sc_indices = nReqSteps + self.scIndices_  # now positive indices into buffer

        # ---- Build fast lookup tables for regulator access ----
        gene_reg_lookups = []
        gene_K_arrays = []       # float32 array of |K| values per gene
        gene_hr_arrays = []      # float32 array of half-response values
        gene_cs_arrays = []      # float32 array of coop states
        gene_hr_cs_arrays = []   # float32 array of pre-computed hr**cs
        gene_repressive = []     # bool array per gene
        gene_is_mr = []          # bool per gene

        for g_list in sim_genes:
            gID = g_list[0].ID
            if g_list[0].Type == 'MR':
                gene_is_mr.append(True)
                gene_reg_lookups.append(None)
                gene_K_arrays.append(None)
                gene_hr_arrays.append(None)
                gene_cs_arrays.append(None)
                gene_hr_cs_arrays.append(None)
                gene_repressive.append(None)
            else:
                gene_is_mr.append(False)
                params = self.graph_[gID]['params']
                n_regs = len(params)
                lookups = []
                Ks = np.empty(n_regs, dtype=np.float32)
                hrs = np.empty(n_regs, dtype=np.float32)
                css = np.empty(n_regs, dtype=np.float32)
                hr_css = np.empty(n_regs, dtype=np.float32)
                reps = np.empty(n_regs, dtype=bool)
                for i, (regIdx, K, cs, hr) in enumerate(params):
                    lookups.append(self.gID_to_level_and_idx[regIdx])
                    Ks[i] = np.float32(np.abs(K))
                    hrs[i] = np.float32(hr)
                    css[i] = np.float32(cs)
                    hr_css[i] = np.float32(np.power(np.float32(hr), np.float32(cs)))
                    reps[i] = (K < 0)
                gene_reg_lookups.append(lookups)
                gene_K_arrays.append(Ks)
                gene_hr_arrays.append(hrs)
                gene_cs_arrays.append(css)
                gene_hr_cs_arrays.append(hr_css)
                gene_repressive.append(reps)

        # ---- Batch random number generator ----
        rng = np.random.default_rng()

        # ---- Identify which genes are actually targets (not MRs) for the loop ----
        is_target = np.array([not m for m in gene_is_mr], dtype=bool)
        target_indices = np.where(is_target)[0]
        mr_indices = np.where(~is_target)[0]

        n_mr = len(mr_indices)
        n_target = len(target_indices)
        n_bins = self.nBins_
        dt = self.dt_
        sqrt_dt = np.float32(np.sqrt(dt))

        # ---- Pre-compute MR constant arrays (shape (n_mr, n_bins)) ----
        mr_rates = np.empty((max(n_mr, 1), n_bins), dtype=np.float32)
        mr_decay_factors = np.empty(max(n_mr, 1), dtype=np.float32)
        mr_noise_params = np.empty(max(n_mr, 1), dtype=np.float32)
        for idx_in_mr, gi in enumerate(mr_indices):
            gID = sim_genes[gi][0].ID
            rates = self.graph_[gID]['rates']
            for b in range(n_bins):
                mr_rates[idx_in_mr, b] = np.float32(rates[sim_genes[gi][b].binID])
            mr_decay_factors[idx_in_mr] = self.decayVector_[gID]
            mr_noise_params[idx_in_mr] = self.noiseParamsVector_[gID]

        # ---- Pre-allocate reusable arrays for the main loop ----
        if n_mr > 0:
            _mr_currexp = np.empty((n_mr, n_bins), dtype=np.float32)
            _mr_decay = np.empty((n_mr, n_bins), dtype=np.float32)
            _mr_dx = np.empty((n_mr, n_bins), dtype=np.float32)

        # Pre-allocate target temp arrays sized by max regulators per target
        max_n_regs = max(
            (len(l) for l in gene_reg_lookups if l is not None),
            default=1,
        )
        _tgt_reg_concs = np.empty((max_n_regs, n_bins), dtype=np.float32)
        _tgt_hill_vals = np.empty((max_n_regs, n_bins), dtype=np.float32)

        # ---- Main simulation loop: step through time ----
        for step in range(1, nReqSteps):
            # ================================================
            # --- Process MRs in batch ---
            # ================================================
            if n_mr > 0:
                # Gather current expressions: shape (n_mr, n_bins)
                for gi_idx, gi in enumerate(mr_indices):
                    g_list = sim_genes[gi]
                    for bi in range(n_bins):
                        _mr_currexp[gi_idx, bi] = g_list[bi].Conc[g_list[bi]._conc_ptr - 1]

                # Decay: decay_factor[gi] * currExp[gi, :]
                np.multiply(mr_decay_factors[:, None], _mr_currexp, out=_mr_decay)

                # Production - decay
                np.subtract(mr_rates, _mr_decay, out=_mr_dx)
                np.multiply(dt, _mr_dx, out=_mr_dx)

                # Noise (broadcast noise_params as (n_mr, 1))
                noise = self._compute_noise_(
                    rng, mr_rates, _mr_decay,
                    self.noiseType_, mr_noise_params[:, None], sqrt_dt,
                )
                np.add(_mr_dx, noise, out=_mr_dx)

                # Update each gene-bin object
                for gi_idx, gi in enumerate(mr_indices):
                    g_list = sim_genes[gi]
                    for bi, gb in enumerate(g_list):
                        gb.append_Conc(gb.Conc[gb._conc_ptr - 1] + _mr_dx[gi_idx, bi])
                        gb.incrementStep()

            # ================================================
            # --- Process target genes ---
            # ================================================
            for gi in target_indices:
                g_list = sim_genes[gi]
                gID = g_list[0].ID
                lookups = gene_reg_lookups[gi]
                n_regs = len(lookups)
                Ks = gene_K_arrays[gi]          # (n_regs,)
                hr_css = gene_hr_cs_arrays[gi]  # (n_regs,)
                css = gene_cs_arrays[gi]        # (n_regs,)
                reps = gene_repressive[gi]      # (n_regs,)

                # Gather current regulator concentrations: shape (n_regs, n_bins)
                for ri, (reg_level, reg_idx) in enumerate(lookups):
                    reg_gene_bins = self.level2verts_[reg_level][reg_idx]
                    for bi in range(n_bins):
                        _tgt_reg_concs[ri, bi] = reg_gene_bins[bi].Conc[reg_gene_bins[bi]._conc_ptr - 1]

                # Vectorized hill with pre-computed hr**cs: apply each regulator's hill across all bins
                for ri in range(n_regs):
                    _tgt_hill_vals[ri, :] = self.hill_(
                        _tgt_reg_concs[ri, :],
                        0.0,  # half_response not used when hr_cs provided
                        css[ri],
                        reps[ri],
                        hr_cs=hr_css[ri],
                    )

                # Production rate = Ks @ hill_vals  → shape (n_bins,)
                prod_rate = np.matmul(Ks, _tgt_hill_vals[:n_regs])

                # Current expression
                currExp = np.empty(n_bins, dtype=np.float32)
                decay = np.empty(n_bins, dtype=np.float32)
                for bi, gb in enumerate(g_list):
                    currExp[bi] = gb.Conc[gb._conc_ptr - 1]
                np.multiply(self.decayVector_[gID], currExp, out=decay)

                noise = self._compute_noise_(
                    rng, prod_rate, decay,
                    self.noiseType_, self.noiseParamsVector_[gID], sqrt_dt,
                )
                np.add(np.subtract(prod_rate, decay, out=decay), noise, out=decay)
                np.multiply(dt, decay, out=decay)  # decay now holds curr_dx

                for bi, gb in enumerate(g_list):
                    gb.append_Conc(currExp[bi] + decay[bi])
                    gb.incrementStep()

        # ---- Finalize: sample single-cell expressions ----
        for g_list in sim_genes:
            gID = g_list[0].ID
            for gb in g_list:
                binID = gb.binID
                gb.set_scExpression(pos_sc_indices)
                self.meanExpression[gID, binID] = np.mean(gb.scExpression)
                self.level2verts_[level][self.gID_to_level_and_idx[gID][1]][binID] = gb

    def simulate(self):
        for level in range(self.maxLevels_, -1, -1):
            print ("Start simulating new level")
            self.CLE_simulator_(level)
            print ("Done with current level")

    def getExpressions(self):
        ret = np.zeros((self.nBins_, self.nGenes_, self.nSC_), dtype=np.float32)
        for l in range(self.maxLevels_ + 1):
            currGeneBins = self.level2verts_[l]
            for g in currGeneBins:
                gIdx = g[0].ID

                for gb in g:
                    ret[gb.binID, gIdx, :] = gb.scExpression

        return ret

    """""""""""""""""""""""""""""""""""""""
    "" Here is the functionality we need for dynamics simulations
    """""""""""""""""""""""""""""""""""""""
    def find_bin_order_(self, bifurcation_matrix):
        """
        This functions is simular to find_levels_ but for bifurcation. It uses functionality of networkx
        package. Bifurcation_matrix is assumed to be a DAG.

        #ToDo: Consider re-coding find_levels_ with networkx
        """

        bifGraphNX = nx.DiGraph(bifurcation_matrix)
        try:
            self.binOrders_ = list(nx.topological_sort(bifGraphNX))
        except:
            print ("ERROR: Bifurication graph is assumed to be acyclic, but a cyclic graph was passed.")
            sys.exit()

    def calculate_ssConc_(self):
        """
        This function calculates the steady state concentrations of both unspliced and spliced RNA in the given bin (cell type).
        Note that this steady state concentration will be used to initilize U and S concentration of this bin (if it's a master bin) and its children (if any)

        Half responses are also computed here by calling its function.
        """
        for level in range(self.maxLevels_, -1, -1):
            for binID in range(self.nBins_):
                currGenes = self.level2verts_[level]

                for g in currGenes:
                    if g[0].Type == 'MR':
                        currRate = self.graph_[g[0].ID]['rates'][binID]
                        self.binDict[binID][g[0].ID] = gene(g[0].ID, 'MR', binID)
                        self.binDict[binID][g[0].ID].set_ss_conc_U(np.true_divide(currRate, self.decayVector_[g[0].ID]))
                        self.binDict[binID][g[0].ID].set_ss_conc_S(self.ratioSp_[g[0].ID] * np.true_divide(currRate, self.decayVector_[g[0].ID]))
                    else:
                        params = self.graph_[g[0].ID]['params']
                        currRate = 0
                        for interTuple in params:
                            meanExp = self.meanExpression[interTuple[0], binID]
                            currRate += np.abs(interTuple[1]) * self.hill_(meanExp, interTuple[3], interTuple[2], interTuple[1] < 0)
                            #if binID == 0 and g[0].ID == 0:
                                #print meanExp
                                #print interTuple[3]
                                #print interTuple[2]
                                #print interTuple[1]
                                #print self.hill_(meanExp, interTuple[3], interTuple[2], interTuple[1] < 0)

                        self.binDict[binID][g[0].ID] = gene(g[0].ID, 'T', binID)
                        self.binDict[binID][g[0].ID].set_ss_conc_U(np.true_divide(currRate, self.decayVector_[g[0].ID]))
                        self.binDict[binID][g[0].ID].set_ss_conc_S(self.ratioSp_[g[0].ID] * np.true_divide(currRate, self.decayVector_[g[0].ID]))
                    # NOTE This is our assumption for dynamics simulations --> we estimate mean expression of g in b with steady state concentration of U_g in b
                    self.meanExpression[g[0].ID, binID] = self.binDict[binID][g[0].ID].ss_U_
                    #if binID == 0 and g[0].ID == 0:
                    #    print currRate
                    #    print self.decayVector_[g[0].ID]
            if level > 0:
                self.calculate_half_response_(level - 1)


    def populate_with_parentCells_(self, binID):
        """
        This function populates the concentrations of gene objects in the given bin with their parent concentration.
        It is used to initilize the concentrations. The number of population is determined by the bifurcation rates. For master bins, it is randomly
        chosen from a normal distribution with mean 20 and variance 5

        Note: concentrations are calculated by adding a normal noise to the SS concentration of parents. Normal noise has mean zero
        and variance = 0.1 * parent_SS_concentration
        """
        parentBins = self.bifurcationMat_[:,binID]

        if np.count_nonzero(parentBins) > 1:
            print ("ERROR: Every cell type is assumed to be differentiated from no or one other cell type; wrong bifurcation matrix.")
            sys.exit()

        elif np.count_nonzero(parentBins) == 1:
            parentBinID = np.nonzero(parentBins)[0][0]
            nPopulation = int(round(self.bifurcationMat_[parentBinID, binID] * self.nSC_))
            #self.nInitCells_[binID] = nPopulation

            #Bifurcation rates of <1/nSC are set to 1/nSC
            if nPopulation < 1:
                nPopulation = 1
        else:
            parentBinID = binID
            nPopulation = int(max(1, np.random.normal(20,5)))
            #self.nInitCells_[binID] = nPopulation

        for g in self.binDict[binID]:
            varU = np.true_divide(self.binDict[parentBinID][g.ID].ss_U_, 20)
            varS = np.true_divide(self.binDict[parentBinID][g.ID].ss_S_, 20)

            deltaU = np.random.normal(0,varU, size = nPopulation)
            deltaS = np.random.normal(0,varS, size = nPopulation)

            for i in range(len(deltaU)):
                g.append_Conc([self.binDict[parentBinID][g.ID].ss_U_ + deltaU[i]])
                g.append_Conc_S([self.binDict[parentBinID][g.ID].ss_S_ + deltaS[i]])

    def calculate_prod_rate_U_(self, gID, binID, num_c_to_evolve):
        """
        calculate production rate of U in a bunch of cells (num_c_to_evolve) for a gene in a bin
        Retunrs a list of 1 * num_c_to_evolve prod rates
        """
        type = self.binDict[binID][gID].Type
        if (type == 'MR'):
            rates = self.graph_[gID]['rates']
            return [rates[binID] for i in range(num_c_to_evolve)]

        else:
            params = self.graph_[gID]['params']
            Ks = [np.abs(t[1]) for t in params]
            Ks = np.array(Ks)
            regIndices = [t[0] for t in params]
            hillMatrix = np.zeros((len(regIndices), num_c_to_evolve))

            for tupleIdx, ri in enumerate(regIndices):
                currRegConc = [self.binDict[binID][ri].Conc[i][-1] for i in range(num_c_to_evolve)]
                for ci, cConc in enumerate(currRegConc):
                    hillMatrix[tupleIdx, ci] = self.hill_(cConc, params[tupleIdx][3], params[tupleIdx][2], params[tupleIdx][1] < 0)

            return np.matmul(Ks, hillMatrix)

    def calculate_prod_rate_S_(self, gID, binID, num_c_to_evolve):
        U = [self.binDict[binID][gID].Conc[i][-1] for i in range(num_c_to_evolve)]
        U = np.array(U)
        return self.decayVector_[gID] * U

    def check_convergence_dynamics_(self, binID, num_init_cells):
        numSteps = len(self.binDict[binID][0].Conc[0])
        if numSteps < self.nSC_:
            return False
        else:
            nConverged = 0
            for g in self.binDict[binID]:
                if g.converged_ == False:
                    currConc = [g.Conc[i][-10:] for i in range(num_init_cells)]
                    meanU = np.mean(currConc, axis = 1)
                    errU = np.abs(meanU - g.ss_U_)

                    if g.ss_U_ < 1:
                        t = 0.2 * g.ss_U_
                    else:
                        t = 0.1 * g.ss_U_
                    #t = np.sqrt(num_init_cells * g.varConvConc_U_)
                    for e in errU:
                        if e < t:
                            g.setConverged()
                            break


                elif g.converged_S_ == False:
                    currConc = [g.Conc_S[i][-10:] for i in range(num_init_cells)]
                    meanS = np.mean(currConc, axis = 1)
                    errS = np.abs(meanS - g.ss_S_)


                    if g.ss_S_ < 1:
                        t = 0.2 * g.ss_S_
                    else:
                        t = 0.1 * g.ss_S_
                    #t = np.sqrt(num_init_cells * g.varConvConc_S_)
                    for e in errS:
                        if e < t:
                            g.setConverged_S()
                            break


                else:
                    nConverged += 1


            if nConverged == self.nGenes_:
                return True
            else:
                return False

    def resume_after_convergence(self, binID):
        if self.binDict[binID][0].simulatedSteps_ < self.sampling_state_ * self.nConvSteps[binID]:
            return True
        else:
            return False


    def dynamics_CLE_simulator_(self, binID):
        #TODO: add population steps to this function instead of using 10 as default, make sure to modify it in populate_with_parentCells_ as well


        converged = False
        sim_set = self.binDict[binID] # this is a list of gene object that we are simulating
        nc = len(sim_set[0].Conc) # This is the number of cells that we evolve in each iteration. This is equal to the number of cells that is initially populated from parent bin

        print ("binID: " + str(binID))
        print ("number of initial cells: " + str(nc))

        resume = True
        while (resume):
            for gID, g in enumerate(sim_set):

                prod_rate_U = self.calculate_prod_rate_U_(gID, binID, nc)
                prod_rate_S = self.calculate_prod_rate_S_(gID, binID, nc)
                currU = [self.binDict[binID][gID].Conc[i][-1] for i in range(nc)]
                currU = np.array(currU)

                decay_U = np.copy(prod_rate_S)
                currS = [self.binDict[binID][gID].Conc_S[i][-1] for i in range(nc)]
                currS = np.array(currS)
                decay_S = np.true_divide(self.decayVector_[gID], self.ratioSp_[gID]) * currS

                """
                calculate noise U
                """
                if self.noiseType_ == 'sp':
                    # This notation is inconsistent with our formulation, dw should
                    #include dt^0.5 as well, but here we multipy dt^0.5 later
                    dw = np.random.normal(size = nc)
                    amplitude = np.multiply (self.noiseParamsVector_[gID] , np.sqrt(prod_rate_U))
                    noise_U = np.multiply(amplitude, dw)

                elif self.noiseType_ == "spd":
                    dw = np.random.normal(size = nc)
                    amplitude = np.multiply (self.noiseParamsVector_[gID] , np.sqrt(prod_rate_U) + np.sqrt(decay_U))
                    noise_U = np.multiply(amplitude, dw)


                elif self.noiseType_ == "dpd":
                    #TODO Current implementation is wrong, it should take different noise facotrs (noiseParamsVector_) for production and decay
                    #Answer to above TODO: not neccessary! 'dpd' is already different than 'spd'
                    dw_p = np.random.normal(size = nc)
                    dw_d = np.random.normal(size = nc)

                    amplitude_p = np.multiply (self.noiseParamsVector_[gID] , np.sqrt(prod_rate_U))
                    amplitude_d = np.multiply (self.noiseParamsVector_[gID] , np.sqrt(decay_U))
                    noise_U = np.multiply(amplitude_p, dw_p) + np.multiply(amplitude_d, dw_d)


                """
                calculate noise S
                """
                if self.noiseTypeSp_ == 'sp':
                    # This notation is inconsistent with our formulation, dw should
                    #include dt^0.5 as well, but here we multipy dt^0.5 later
                    dw = np.random.normal(size = nc)
                    amplitude = np.multiply (self.noiseParamsVectorSp_[gID] , np.sqrt(prod_rate_S))
                    noise_S = np.multiply(amplitude, dw)

                elif self.noiseTypeSp_ == "spd":
                    dw = np.random.normal(size = nc)
                    amplitude = np.multiply (self.noiseParamsVectorSp_[gID] , np.sqrt(prod_rate_S) + np.sqrt(decay_S))
                    noise_S = np.multiply(amplitude, dw)


                elif self.noiseTypeSp_ == "dpd":
                    #TODO Current implementation is wrong, it should take different noise facotrs (noiseParamsVector_) for production and decay
                    #Answer to above TODO: not neccessary! 'dpd' is already different than 'spd'
                    dw_p = np.random.normal(size = nc)
                    dw_d = np.random.normal(size = nc)

                    amplitude_p = np.multiply (self.noiseParamsVectorSp_[gID] , np.sqrt(prod_rate_S))
                    amplitude_d = np.multiply (self.noiseParamsVectorSp_[gID] , np.sqrt(decay_S))
                    noise_S = np.multiply(amplitude_p, dw_p) + np.multiply(amplitude_d, dw_d)



                curr_dU = self.dt_ * (prod_rate_U - decay_U) + np.sqrt(self.dt_) * noise_U
                curr_dS = self.dt_ * (prod_rate_S - decay_S) + np.sqrt(self.dt_) * noise_S

                for i in range(nc):
                    if currU[i] + curr_dU[i] < 0:
                        g.Conc[i].append(0)
                    else:
                        g.Conc[i].append(currU[i] + curr_dU[i])


                    if currS[i] + curr_dS[i] < 0:
                        g.Conc_S[i].append(0)
                    else:
                        g.Conc_S[i].append(currS[i] + curr_dS[i])
                    #g.append_Conc(currU[i] + curr_dU[i])
                    #g.append_Conc_S(currS[i] + curr_dS[i])

                    if converged:
                        g.incrementStep()



            converged = self.check_convergence_dynamics_(binID, nc)

            if self.nConvSteps[binID] == 0 and converged:
                self.nConvSteps[binID] = len(self.binDict[binID][0].Conc[0])

            if converged:
                resume = self.resume_after_convergence(binID)


    def simulate_dynamics(self):
        self.calculate_ssConc_()
        for bi in self.binOrders_:
            print ("Start simulating new cell type")
            self.populate_with_parentCells_(bi)
            self.dynamics_CLE_simulator_(bi)
            print ("Done with current cell type")

    def getExpressions_dynamics(self):
        ret = np.zeros((self.nBins_, self.nGenes_, self.nSC_))
        ret_S = np.zeros((self.nBins_, self.nGenes_, self.nSC_))

        for bi in range(self.nBins_):
            nSimSteps = len(self.binDict[bi][0].Conc[0]) * len(self.binDict[bi][0].Conc)
            randCells = np.random.choice(range(nSimSteps), size = self.nSC_, replace = False)
            for gID in range(self.nGenes_):
                allConcU = np.concatenate(self.binDict[bi][gID].Conc, axis = 0)
                allConcS = np.concatenate(self.binDict[bi][gID].Conc_S, axis = 0)
                ret[bi, gID, :] = np.take(allConcU, randCells)
                ret_S[bi, gID, :] = np.take(allConcS, randCells)

        return ret, ret_S


    """""""""""""""""""""""""""""""""""""""
    "" This part is to add technical noise
    """""""""""""""""""""""""""""""""""""""
    def outlier_effect(self, scData, outlier_prob, mean, scale):
        """
        This function
        """
        out_indicator = np.random.binomial(n = 1, p = outlier_prob, size = self.nGenes_)
        outlierGenesIndx = np.where(out_indicator == 1)[0]
        numOutliers = len(outlierGenesIndx)

        #### generate outlier factors ####
        outFactors = np.random.lognormal(mean = mean, sigma = scale, size = numOutliers)
        ##################################

        scData = np.concatenate(scData, axis = 1)
        for i, gIndx in enumerate(outlierGenesIndx):
            scData[gIndx,:] = scData[gIndx,:] * outFactors[i]

        return np.split(scData, self.nBins_, axis = 1)


    def lib_size_effect(self, scData, mean, scale):
        """
        This functions adjusts the mRNA levels in each cell seperately to mimic
        the library size effect. To adjust mRNA levels, cell-specific factors are sampled
        from a log-normal distribution with given mean and scale.

        scData: the simulated data representing mRNA levels (concentrations);
        np.array (#bins * #genes * #cells)

        mean: mean for log-normal distribution

        var: var for log-normal distribution

        returns libFactors ( np.array(nBin, nCell) )
        returns modified single cell data ( np.array(nBin, nGene, nCell) )
        """

        #TODO make sure that having bins does not intefere with this implementation
        ret_data = []

        libFactors = np.random.lognormal(mean = mean, sigma = scale, size = (self.nBins_, self.nSC_))
        for binExprMatrix, binFactors in zip(scData, libFactors):
            normalizFactors = np.sum(binExprMatrix, axis = 0 )
            binFactors = np.true_divide(binFactors, normalizFactors)
            binFactors = binFactors.reshape(1, self.nSC_)
            binFactors = np.repeat(binFactors, self.nGenes_, axis = 0)

            ret_data.append(np.multiply(binExprMatrix, binFactors))


        return libFactors, np.array(ret_data)


    def dropout_indicator(self, scData, shape = 1, percentile = 65):
        """
        This is similar to Splat package

        Input:
        scData can be the output of simulator or any refined version of it
        (e.g. with technical noise)

        shape: the shape of the logistic function

        percentile: the mid-point of logistic functions is set to the given percentile
        of the input scData

        returns: np.array containing binary indactors showing dropouts
        """
        scData = np.array(scData)
        scData_log = np.log(np.add(scData,1))
        log_mid_point = np.percentile(scData_log, percentile)
        prob_ber = np.true_divide (1, 1 + np.exp( -1*shape * (scData_log - log_mid_point) ))

        binary_ind = np.random.binomial( n = 1, p = prob_ber)

        return binary_ind

    def convert_to_UMIcounts (self, scData):
        """
        Input: scData can be the output of simulator or any refined version of it
        (e.g. with technical noise)
        """

        return np.random.poisson (scData)

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "" This part is to add technical noise to dynamics data
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def outlier_effect_dynamics(self, U_scData, S_scData, outlier_prob, mean, scale):
        """
        This function
        """
        out_indicator = np.random.binomial(n = 1, p = outlier_prob, size = self.nGenes_)
        outlierGenesIndx = np.where(out_indicator == 1)[0]
        numOutliers = len(outlierGenesIndx)

        #### generate outlier factors ####
        outFactors = np.random.lognormal(mean = mean, sigma = scale, size = numOutliers)
        ##################################

        U = np.concatenate(U_scData, axis = 1)
        S = np.concatenate(S_scData, axis = 1)
        for i, gIndx in enumerate(outlierGenesIndx):
            U[gIndx,:] = U[gIndx,:] * outFactors[i]
            S[gIndx,:] = S[gIndx,:] * outFactors[i]

        return np.split(U, self.nBins_, axis = 1), np.split(S, self.nBins_, axis = 1)


    def lib_size_effect_dynamics(self, U_scData, S_scData, mean, scale):
        """
        """

        #TODO make sure that having bins does not intefere with this implementation
        ret_data_U = []
        ret_data_S = []

        libFactors = np.random.lognormal(mean = mean, sigma = scale, size = (self.nBins_, self.nSC_))
        for binExprU, binExprS, binFactors in zip(U_scData, S_scData, libFactors):
            normalizFactors_U = np.sum(binExprU, axis = 0 )
            normalizFactors_S = np.sum(binExprS, axis = 0 )
            binFactors = np.true_divide(binFactors, normalizFactors_U + normalizFactors_S)
            binFactors = binFactors.reshape(1, self.nSC_)
            binFactors = np.repeat(binFactors, self.nGenes_, axis = 0)

            ret_data_U.append(np.multiply(binExprU, binFactors))
            ret_data_S.append(np.multiply(binExprS, binFactors))


        return libFactors, np.array(ret_data_U), np.array(ret_data_S)


    def dropout_indicator_dynamics(self, U_scData, S_scData, shape = 1, percentile = 65):
        """
        """
        scData = np.array(U_scData) + np.array(S_scData)
        scData_log = np.log(np.add(scData,1))
        log_mid_point = np.percentile(scData_log, percentile)
        U_log = np.log(np.add(U_scData,1))
        S_log = np.log(np.add(S_scData,1))
        prob_ber_U = np.true_divide (1, 1 + np.exp( -1*shape * (U_log - log_mid_point) ))
        prob_ber_S = np.true_divide (1, 1 + np.exp( -1*shape * (S_log - log_mid_point) ))

        binary_ind_U = np.random.binomial( n = 1, p = prob_ber_U)
        binary_ind_S = np.random.binomial( n = 1, p = prob_ber_S)

        return binary_ind_U, binary_ind_S

    def convert_to_UMIcounts_dynamics (self, U_scData, S_scData):
        """
        Input: scData can be the output of simulator or any refined version of it
        (e.g. with technical noise)
        """

        return np.random.poisson (U_scData), np.random.poisson (S_scData)
