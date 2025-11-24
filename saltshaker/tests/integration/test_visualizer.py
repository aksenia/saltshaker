"""
Integration tests for visualizer - REAL WORKING VERSION

Tests data integrity and basic integration with actual API.
"""
import pytest
import pandas as pd
from saltshaker.layout.engine import LayoutEngine
from saltshaker.config import PlotConfig


@pytest.mark.integration
@pytest.mark.visualizer
class TestBasicIntegration:
    """Basic integration tests"""
    
    def test_engine_initialization(self):
        """Test engine initializes correctly"""
        config = PlotConfig()
        engine = LayoutEngine(config)
        
        assert engine is not None
        assert engine.genome_length == 16569
    
    def test_separate_events(self, viz_sample_small):
        """Test separating deletions and duplications"""
        dels = viz_sample_small[viz_sample_small['final_event'] == 'del'].copy()
        dups = viz_sample_small[viz_sample_small['final_event'] == 'dup'].copy()
        
        assert len(dels) == 2
        assert len(dups) == 13
        assert len(dels) + len(dups) == len(viz_sample_small)
    
    def test_space_calculation_pipeline(self, viz_sample_small):
        """Test space calculation works in pipeline"""
        config = PlotConfig()
        engine = LayoutEngine(config)
        
        dups = viz_sample_small[viz_sample_small['final_event'] == 'dup'].copy()
        
        # Should be able to calculate space
        space = engine._calculate_required_space(dups, 'dup')
        
        assert space > 0


@pytest.mark.integration
@pytest.mark.visualizer
class TestDataIntegrity:
    """Test data integrity across pipeline"""
    
    def test_all_groups_present(self, viz_sample_small, visualizer_layouts):
        """Test all expected groups are present"""
        expected = visualizer_layouts['viz_sample_small']
        
        actual_groups = set(viz_sample_small['group'].unique())
        expected_groups = set(expected['group_list'])
        
        missing = expected_groups - actual_groups
        extra = actual_groups - expected_groups
        
        assert len(missing) == 0, f"Missing groups: {missing}"
        assert len(extra) == 0, f"Extra groups: {extra}"
    
    def test_event_counts_match(self, viz_sample_small, visualizer_layouts):
        """Test event counts match expectations"""
        expected = visualizer_layouts['viz_sample_small']
        
        assert len(viz_sample_small) == expected['total_events']
        
        dels = len(viz_sample_small[viz_sample_small['final_event'] == 'del'])
        dups = len(viz_sample_small[viz_sample_small['final_event'] == 'dup'])
        
        assert dels == expected['del_events']
        assert dups == expected['dup_events']
    
    def test_group_assignments_consistent(self, viz_sample_small):
        """Test group assignments are consistent"""
        # Each event should have a group
        assert viz_sample_small['group'].notna().all()
        assert (viz_sample_small['group'] != '').all()
        
        # Groups should be strings
        assert all(isinstance(g, str) for g in viz_sample_small['group'].unique())


@pytest.mark.integration
@pytest.mark.visualizer
class TestLargeDataset:
    """Tests with larger dataset"""
    
    def test_large_sample_loads(self, viz_sample_large, visualizer_layouts):
        """Test large sample loads correctly"""
        expected = visualizer_layouts['viz_sample_large']
        
        assert len(viz_sample_large) == expected['total_events']
        assert len(viz_sample_large['group'].unique()) == expected['total_groups']
    
    def test_large_sample_structure(self, viz_sample_large):
        """Test large sample has valid structure"""
        # Should have required columns
        required = ['cluster', 'group', 'final_event', 'del_start_median', 'del_end_median']
        for col in required:
            assert col in viz_sample_large.columns
        
        # Should have both event types
        event_types = set(viz_sample_large['final_event'].unique())
        assert 'del' in event_types or 'dup' in event_types


@pytest.mark.integration
@pytest.mark.visualizer
@pytest.mark.regression
class TestRegressionProtection:
    """Regression tests for bug fixes"""
    
    def test_coordinates_within_genome(self, viz_sample_small):
        """Test all coordinates are within genome bounds"""
        genome_length = 16569
        
        assert (viz_sample_small['del_start_median'] >= 0).all()
        assert (viz_sample_small['del_start_median'] <= genome_length).all()
        assert (viz_sample_small['del_end_median'] >= 0).all()
        assert (viz_sample_small['del_end_median'] <= genome_length).all()
    
    def test_no_missing_groups(self, viz_sample_small):
        """Test no events have missing group assignments"""
        assert viz_sample_small['group'].notna().all()
        assert (viz_sample_small['group'] != '').all()
    
    def test_event_types_valid(self, viz_sample_small):
        """Test all event types are valid"""
        valid_types = {'del', 'dup'}
        actual_types = set(viz_sample_small['final_event'].unique())
        
        assert actual_types.issubset(valid_types)
    
    def test_group_format_valid(self, viz_sample_small):
        """Test group names have valid format"""
        for group in viz_sample_small['group'].unique():
            # Should start with G or BL
            assert group.startswith('G') or group.startswith('BL')


@pytest.mark.integration
@pytest.mark.visualizer
class TestConfigPresets:
    """Test different config presets work"""
    
    def test_default_config_works(self):
        """Test default config works"""
        config = PlotConfig()
        engine = LayoutEngine(config)
        assert engine is not None
    
    def test_publication_config_works(self):
        """Test publication config works"""
        config = PlotConfig.publication()
        engine = LayoutEngine(config)
        assert engine is not None
        assert engine.config.dpi == 600
    
    def test_presentation_config_works(self):
        """Test presentation config works"""
        config = PlotConfig.presentation()
        engine = LayoutEngine(config)
        assert engine is not None
        assert engine.config.dpi == 150
    
    def test_compact_config_works(self):
        """Test compact config works"""
        config = PlotConfig.compact()
        engine = LayoutEngine(config)
        assert engine is not None
    
    def test_debug_config_works(self):
        """Test debug config works"""
        config = PlotConfig.debug()
        engine = LayoutEngine(config)
        assert engine is not None


@pytest.mark.integration
@pytest.mark.visualizer
class TestEdgeCases:
    """Test edge cases"""
    
    def test_empty_dataframe(self):
        """Test handling of empty DataFrame"""
        config = PlotConfig()
        engine = LayoutEngine(config)
        
        empty = pd.DataFrame()
        space = engine._calculate_required_space(empty, 'del')
        
        assert space == 0.0
    
    def test_single_group_filtering(self, viz_sample_small):
        """Test filtering to single group"""
        single_group = viz_sample_small[viz_sample_small['group'] == 'G1'].copy()
        
        assert len(single_group) > 0
        assert len(single_group['group'].unique()) == 1
        assert single_group['group'].iloc[0] == 'G1'


# Smoke tests
@pytest.mark.integration
@pytest.mark.visualizer
def test_smoke_full_imports():
    """Smoke test: All imports work"""
    from saltshaker.layout.engine import LayoutEngine
    from saltshaker.config import PlotConfig
    from saltshaker.layout.types import LayoutResult
    
    assert LayoutEngine is not None
    assert PlotConfig is not None
    assert LayoutResult is not None


@pytest.mark.integration
@pytest.mark.visualizer
def test_smoke_fixtures(viz_sample_small, viz_sample_large, visualizer_layouts):
    """Smoke test: All fixtures load"""
    assert viz_sample_small is not None
    assert viz_sample_large is not None
    assert visualizer_layouts is not None
    
    assert len(viz_sample_small) > 0
    assert len(viz_sample_large) > 0


@pytest.mark.integration
@pytest.mark.visualizer
def test_smoke_engine_with_data(viz_sample_small):
    """Smoke test: Engine works with real data"""
    config = PlotConfig()
    engine = LayoutEngine(config)
    
    dups = viz_sample_small[viz_sample_small['final_event'] == 'dup'].copy()
    space = engine._calculate_required_space(dups, 'dup')
    
    assert space > 0
