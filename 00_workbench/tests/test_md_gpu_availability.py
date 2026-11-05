import pytest
from research import md_queue


def test_md_waits_for_transient_gpu_occupancy(monkeypatch,tmp_path):
    snapshots=iter([{'gpus':[{'allocatable':False}]},
                    {'gpus':[{'allocatable':True,'uuid':'public-test-gpu'}]}])
    monkeypatch.setattr(md_queue,'detect',lambda _:next(snapshots))
    monkeypatch.setattr(md_queue.time,'sleep',lambda _:None)
    assert md_queue.available_gpus(tmp_path,1)[0]['uuid']=='public-test-gpu'


def test_busy_gpu_is_not_queued_or_preempted(monkeypatch,tmp_path):
    clock=iter([0,11])
    monkeypatch.setattr(md_queue.time,'monotonic',lambda:next(clock))
    monkeypatch.setattr(md_queue,'detect',lambda _:{'gpus':[],
        'scheduler_gpu_allocation':{'CUDA_VISIBLE_DEVICES':''}})
    with pytest.raises(ValueError,match='task was not queued'):
        md_queue.available_gpus(tmp_path,1)
