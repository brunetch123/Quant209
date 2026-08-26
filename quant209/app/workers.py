"""Background work: the pipeline and MSConvert, off the UI thread, reported by signal.

The rule this module exists to enforce
---------------------------------------
**Nothing here ever touches a widget.** A ``QRunnable`` runs on a ``QThreadPool`` thread,
and Qt widgets may only be read or written from the thread that created them. Doing it
anyway does not raise. It corrupts, usually much later, in a repaint, as a crash with no
relation to the code that caused it. The workers' only channel back to the UI is therefore a
signal, and the UI's only channel in is a cancellation flag.

Why the signals live on a separate object
------------------------------------------
``QRunnable`` is not a ``QObject`` and cannot declare signals. The standard remedy is a
small ``QObject`` holder, constructed on the UI thread (in ``__init__``, which the UI thread
runs) and therefore *owned* by it. Emitting one of its signals from the pool thread is then
a cross-thread emit, which Qt queues and delivers on the UI thread's event loop. It does
so automatically, because ``Qt.AutoConnection`` compares the emitting thread with the
receiver's. That is why ``signals`` must not be created lazily inside :meth:`run`.

``setAutoDelete(False)``
-------------------------
Qt deletes an auto-delete runnable's C++ object as soon as ``run`` returns, which leaves the
Python wrapper pointing at freed memory. The caller *is* still holding it, to read
``worker.signals`` or call ``worker.cancel()``. Python's own reference counting is the right
lifetime owner here, so auto-delete is switched off in ``__init__``.

Cancellation
-------------
A run is cancelled by setting a flag. The flag is checked inside the progress callback,
which the pipeline calls at every checkpoint, once per file during integration and once per
step elsewhere. On the first check after the flag is set, the callback raises
:class:`_Cancelled`. The pipeline then unwinds through its own ``finally`` blocks, releasing
chromatograms properly, and the partial result is abandoned rather than reported. The user
sees the run stop at the next file boundary rather than mid-integration. Interrupting an
integration would leave a file half-read for no benefit, since a cancelled run has no
partial answer worth keeping.

A cancelled run emits :attr:`WorkerSignals.cancelled`. It never emits ``error``, because
the user asked for this and it is not a failure, and never ``finished``, because there is no
result.
"""

from __future__ import annotations

import threading
import traceback
from pathlib import Path
from typing import Any

from PySide6.QtCore import QObject, QRunnable, Signal

#: Phase name the conversion worker reports under, so a single progress widget can render
#: both workers. ``core.pipeline``'s own phases are "calibrate", "integrate", "quantify".
PHASE_CONVERT = "convert"


class _Cancelled(Exception):
    """Raised inside the progress callback to unwind a cancelled run. Never escapes."""


class CancelToken:
    """A thread-safe "stop at the next checkpoint" flag, shared UI thread -> worker thread.

    An ``Event`` rather than a bare ``bool`` so the write is unambiguously visible to the
    worker thread, and so a caller can hold the token independently of the worker (the
    Process page keeps one across a retry).
    """

    def __init__(self) -> None:
        self._event = threading.Event()

    def cancel(self) -> None:
        self._event.set()

    @property
    def cancelled(self) -> bool:
        return self._event.is_set()

    def raise_if_cancelled(self) -> None:
        if self._event.is_set():
            raise _Cancelled()

    def __bool__(self) -> bool:  # pragma: no cover. Convenience for `if token:`
        return self._event.is_set()


class WorkerSignals(QObject):
    """The one channel from a worker thread back to the UI.

    ``progress`` carries ``core.pipeline``'s own callback contract unchanged, which is
    ``(phase, i, n, label)`` with ``i`` 0-based. The GUI and the CLI therefore render the
    same events and neither has to translate.
    """

    progress = Signal(str, int, int, str)
    #: The worker's return value: a ``RunResult`` or a ``list[ConversionResult]``. Typed as
    #: ``object`` so this module never imports either (and so never imports pandas).
    finished = Signal(object)
    #: ``(message, traceback)``, a one-line explanation for the user and the full
    #: traceback for the details expander and the log.
    error = Signal(str, str)
    cancelled = Signal()
    #: Emitted last, whatever the outcome, so a page can re-enable its buttons in one place.
    done = Signal()


class _BaseWorker(QRunnable):
    """Shared plumbing: signal holder, cancel token, and the emit-exactly-one-outcome rule."""

    def __init__(self, cancel_token: CancelToken | None = None) -> None:
        super().__init__()
        self.setAutoDelete(False)  # see the module docstring
        self.signals = WorkerSignals()
        self.token = cancel_token or CancelToken()

    def cancel(self) -> None:
        """Ask the run to stop at its next checkpoint. Safe from the UI thread, any time."""
        self.token.cancel()

    def _work(self) -> Any:  # pragma: no cover. Overridden by each worker.
        raise NotImplementedError

    def run(self) -> None:
        """Runs on a pool thread. Emits exactly one of finished/cancelled/error, then done."""
        try:
            result = self._work()
        except _Cancelled:
            self.signals.cancelled.emit()
        except BaseException as exc:  # noqa: BLE001. A worker thread must never let one escape
            # An exception that reaches the top of a QRunnable.run is printed to stderr and
            # then forgotten. On a desktop application that is a silent failure. Every one of
            # them becomes a message on screen instead.
            self.signals.error.emit(_message(exc), traceback.format_exc())
        else:
            if self.token.cancelled:
                # Cancelled between the last checkpoint and the return: honour the request
                # rather than delivering a result the user has already navigated away from.
                self.signals.cancelled.emit()
            else:
                self.signals.finished.emit(result)
        finally:
            self.signals.done.emit()

    def _progress(self, phase: str, i: int, n: int, label: str) -> None:
        """The pipeline's progress callback: the cancellation checkpoint, and the emit."""
        self.token.raise_if_cancelled()
        self.signals.progress.emit(phase, int(i), int(n), str(label))


def _message(exc: BaseException) -> str:
    """A one-line, analyst-readable rendering of an exception.

    ``str(exc)`` for anything that carries its own explanation (``core.pipeline`` raises
    ``ValueError`` with full sentences), falling back to the class name for the ones that
    do not. An empty message box is worse than a bare ``KeyError``.
    """
    text = str(exc).strip()
    return text or exc.__class__.__name__


class PipelineWorker(_BaseWorker):
    """Runs :func:`quant209.core.pipeline.run_analysis` on a pool thread.

    ``run_analysis`` is imported inside :meth:`_work`, not at module scope, for two
    reasons. It pulls pandas, about 400 ms, and the cold-start budget forbids that before
    the window paints. A test can also monkeypatch ``quant209.core.pipeline.run_analysis``
    and be sure this worker picks up the stub, because the name is resolved at call time.
    Patching ``quant209.app.workers.run_analysis`` would find nothing, since there is no such
    name. ``quant209.cli`` follows the same convention.
    """

    def __init__(
        self,
        entries: list,
        converted_dir: Path | str,
        params: Any = None,
        qaqc: Any = None,
        analytes: list | None = None,
        reference: str | None = None,
        cancel_token: CancelToken | None = None,
    ) -> None:
        super().__init__(cancel_token)
        self.entries = list(entries)
        self.converted_dir = Path(converted_dir)
        self.params = params
        self.qaqc = qaqc
        self.analytes = analytes
        self.reference = reference

    def _work(self) -> Any:
        from quant209.core.pipeline import run_analysis

        kwargs: dict[str, Any] = {"params": self.params, "progress": self._progress}
        # Passed only when set, so this worker never has to know ``run_analysis``'s
        # defaults (the reference policy in particular is a scientific choice that lives
        # in the pipeline, not in the GUI).
        if self.qaqc is not None:
            kwargs["qaqc"] = self.qaqc
        if self.analytes is not None:
            kwargs["analytes"] = self.analytes
        if self.reference is not None:
            kwargs["reference"] = self.reference
        return run_analysis(self.entries, self.converted_dir, **kwargs)


class ReportWorker(_BaseWorker):
    """Writes the analyst-format Excel workbook on a pool thread.

    Off the UI thread for the same reason the pipeline is: ``report.excel.write_report``
    builds one sheet per reported file through openpyxl and takes several seconds on a
    thirty-file run, and a window that stops repainting for several seconds is a window the
    operating system offers to close for you.

    ``write_report`` is imported inside :meth:`_work`, because openpyxl costs about 200 ms
    and is the Export page's whole reason for loading lazily. It is resolved from its own
    module at call time, so a test monkeypatching ``quant209.report.excel.write_report`` is
    picked up.
    Same convention as :class:`PipelineWorker`.
    """

    def __init__(
        self,
        result: Any,
        out_path: Path | str,
        run_name: str | None = None,
        cancel_token: CancelToken | None = None,
    ) -> None:
        super().__init__(cancel_token)
        self.result = result
        self.out_path = Path(out_path)
        self.run_name = run_name

    def _work(self) -> Any:
        from quant209.report.excel import write_report

        return write_report(self.result, self.out_path, run_name=self.run_name)


class ValidationWorker(_BaseWorker):
    """Reads the analyst's ground truth and writes the validation report, on a pool thread.

    Both halves belong to one worker because neither is useful alone: the report cannot be
    written without the truth, and reading twenty analyst workbooks (openpyxl, ~10 s) only
    to discover there is nothing to compare against is the failure the page's disabled state
    already prevents. Raises ``ValueError`` when the folder holds no usable workbook, which
    the page renders as an error banner. The enablement check is a filename scan and can be
    fooled by a workbook that turns out to hold no quantification sheets.
    """

    def __init__(
        self,
        result: Any,
        run_folder: Path | str,
        out_dir: Path | str,
        run_name: str | None = None,
        cancel_token: CancelToken | None = None,
    ) -> None:
        super().__init__(cancel_token)
        self.result = result
        self.run_folder = Path(run_folder)
        self.out_dir = Path(out_dir)
        self.run_name = run_name or self.run_folder.name

    def _work(self) -> Any:
        from quant209.report.truth import read_run_truth
        from quant209.report.validation import generate_validation_report

        truth = read_run_truth(self.run_folder)
        if truth.areas.empty:
            raise ValueError(
                f"no analyst quantification workbook in {self.run_folder} holds "
                "comparable areas, so there is nothing to validate against"
            )
        report = generate_validation_report(self.result, truth, self.out_dir, self.run_name)
        # The directory, not the report object: the page's job is to point at what was
        # written, and every one of the eight files went here.
        return (self.out_dir, report)


class ConvertWorker(_BaseWorker):
    """Runs :func:`quant209.convert.msconvert.convert_many` on a pool thread.

    ``convert_many``'s callback is ``(i, n, name)``, with no phase, because it has only one.
    It is adapted here to the four-argument :attr:`WorkerSignals.progress` contract with
    the phase filled in as :data:`PHASE_CONVERT`, so the Convert and Process pages can
    share one progress widget instead of growing two.

    Cancellation stops the run *between* files: msconvert is a subprocess and killing it
    mid-write is what produces the truncated mzML this project already has to detect
    (``convert.msconvert.is_truncated``). Files converted before the cancellation are
    finished files and stay on disk.
    """

    def __init__(
        self,
        sources: list[Path | str],
        msconvert_exe: Path | str,
        out_dir: Path | str,
        extra_args: tuple[str, ...] = (),
        cancel_token: CancelToken | None = None,
    ) -> None:
        super().__init__(cancel_token)
        self.sources = list(sources)
        self.msconvert_exe = Path(msconvert_exe)
        self.out_dir = Path(out_dir)
        self.extra_args = tuple(extra_args)

    def _work(self) -> Any:
        from quant209.convert.msconvert import convert_many

        def progress(i: int, n: int, name: str) -> None:
            self._progress(PHASE_CONVERT, i, n, name)

        return convert_many(
            self.sources,
            self.msconvert_exe,
            self.out_dir,
            self.extra_args,
            progress=progress,
        )
