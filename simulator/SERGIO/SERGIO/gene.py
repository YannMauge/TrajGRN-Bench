import numpy as np


class gene(object):

    def __init__(self, geneID, geneType, binID = -1):

        """
        geneType: 'MR' master regulator or 'T' target
        bindID is optional
        """

        self.ID = geneID
        self.Type = geneType
        self.binID = binID
        self.Conc = []
        self.Conc_S = []
        self.dConc = []
        self.k = [] #For dynamics simulation it stores k1 to k4 for Rung-Kutta method, list of size 4 * num_c_to_evolve
        self.k_S = [] #For dynamics simulation it stores k1 to k4 for Rung-Kutta method, list of size 4 * num_c_to_evolve
        self.simulatedSteps_ = 0
        self.converged_ = False
        self.converged_S_ = False
        self.ss_U_ = 0 #This is the steady state concentration of Unspliced mRNA
        self.ss_S_ = 0 #This is the steady state concentration of Spliced mRNA
        self._use_buffer = False   # True when pre-allocated array mode is active (steady-state)
        self._conc_ptr = 0         # write cursor into the pre-allocated buffer

    # ------------------------------------------------------------------
    # Pre-allocated buffer API (used by steady-state CLE_simulator_)
    # ------------------------------------------------------------------
    def init_conc_buffer(self, n_steps: int, max_buffer_mb: float = 4096.0):
        """Pre-allocate a float32 array for concentration history.

        After this call, ``append_Conc`` writes into the buffer instead of
        appending to a Python list, and ``Conc[-1]`` / ``Conc[i]`` work as
        before via ``__getitem__`` delegation.

        Parameters
        ----------
        n_steps : int
            Number of simulation steps to allocate.
        max_buffer_mb : float
            Per-object memory limit in MiB.  Raises ``MemoryError`` if the
            requested buffer exceeds it.
        """
        req_mb = n_steps * 4 / (1024 * 1024)  # float32 = 4 bytes
        if req_mb > max_buffer_mb:
            raise MemoryError(
                f"Gene {self.ID} bin {self.binID}: requested buffer "
                f"({req_mb:.1f} MiB) exceeds per-object limit "
                f"({max_buffer_mb:.1f} MiB). "
                f"Reduce sampling_state, n_cells, or increase max_buffer_mb."
            )
        self._use_buffer = True
        self._conc_ptr = 0
        self.Conc = np.empty(n_steps, dtype=np.float32)

    def _buffer_append(self, val):
        if self._conc_ptr >= len(self.Conc):
            raise MemoryError(
                f"Gene {self.ID} bin {self.binID}: concentration buffer "
                f"overflow ({self._conc_ptr} >= {len(self.Conc)}). "
                f"This indicates a bug in step-count estimation."
            )
        self.Conc[self._conc_ptr] = max(val, 0.0)
        self._conc_ptr += 1

    def _buffer_getitem(self, idx):
        """Delegate indexing so that Conc[-1] etc. work on the filled portion."""
        if self._conc_ptr == 0:
            raise IndexError("buffer is empty")
        return self.Conc[:self._conc_ptr][idx]

    def append_Conc(self, currConc):
        if self._use_buffer:
            self._buffer_append(currConc)
            return
        if isinstance(currConc, list):
            if currConc[0] < 0:
                self.Conc.append([0])
            else:
                self.Conc.append(currConc)
        else:
            if currConc < 0:
                self.Conc.append(0)
            else:
                self.Conc.append(currConc)

    def append_Conc_S(self, currConc):
        if isinstance(currConc, list):
            if currConc[0] < 0:
                self.Conc_S.append([0])
            else:
                self.Conc_S.append(currConc)
        else:
            if currConc < 0:
                self.Conc_S.append(0)
            else:
                self.Conc_S.append(currConc)

    def append_dConc (self, currdConc):
        self.dConc.append(currdConc)

    def append_k (self, list_currK):
        self.k.append(list_currK)

    def append_k_S (self, list_currK):
        self.k_S.append(list_currK)

    def del_lastK_Conc(self, K):
        for k in range(K):
            self.Conc.pop(-1)

    def del_lastK_Conc_S(self, K):
        for k in range(K):
            self.Conc_S.pop(-1)

    def clear_Conc (self):
        """
        This method clears all the concentrations except the last one that may
        serve as intial condition for rest of the simulations
        """
        self.Conc = self.Conc[-1:]

    def clear_dConc (self):
        self.dConc = []

    def incrementStep (self):
        self.simulatedSteps_ += 1

    def setConverged (self):
        self.converged_ = True

    def setConverged_S (self):
        self.converged_S_ = True

    def set_scExpression(self, list_indices):
        """
        selects input indices from self.Conc and form sc Expression
        """
        if self._use_buffer:
            filled = self.Conc[:self._conc_ptr]
            self.scExpression = filled[list_indices]
        else:
            self.scExpression = np.array(self.Conc)[list_indices]

    def set_ss_conc_U(self, u_ss):
        if u_ss < 0:
            u_ss = 0

        self.ss_U_ = u_ss

    def set_ss_conc_S(self, s_ss):
        if s_ss < 0:
            s_ss = 0

        self.ss_S_ = s_ss

    def clear_k(self):
        self.k = []

    def clear_k_S(self):
        self.k_S = []
